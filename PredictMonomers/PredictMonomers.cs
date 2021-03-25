using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Warp;
using Warp.Tools;


namespace NearestNeighbor
{
    struct Particle
    {
        public float3 pos;
        public Matrix3 orientation;
        public int idx;
        public String micrographName;
        public int micrographIdx;
        public Particle(float3 p, float3 o, int i, String mName, int mIdx)
        {
            idx = i;
            pos = p;
            orientation = Matrix3.Euler(o);
            micrographName = mName;
            micrographIdx = mIdx;
        }
    }




    class NearestNeigbour
    {
        static float MatrixDist(Matrix3 mat1, Matrix3 mat2)
        {
            return (float)Math.Sqrt(Math.Pow(mat1.M11 - mat2.M11, 2) + Math.Pow(mat1.M12 - mat2.M12, 2) + Math.Pow(mat1.M13 - mat2.M13, 2)
                + Math.Pow(mat1.M21 - mat2.M21, 2) + Math.Pow(mat1.M22 - mat2.M22, 2) + Math.Pow(mat1.M23 - mat2.M23, 2)
                + Math.Pow(mat1.M31 - mat2.M31, 2) + Math.Pow(mat1.M32 - mat2.M32, 2) + Math.Pow(mat1.M33 - mat2.M33, 2));
        }

        static float3[] GetRelionAngstromOffsets(Star starFile)
        {
            float[] X = starFile.HasColumn("rlnOriginXAngst") ? starFile.GetColumn("rlnOriginXAngst").Select(v => float.Parse(v, CultureInfo.InvariantCulture)).ToArray() : new float[starFile.RowCount];
            float[] Y = starFile.HasColumn("rlnOriginYAngst") ? starFile.GetColumn("rlnOriginYAngst").Select(v => float.Parse(v, CultureInfo.InvariantCulture)).ToArray() : new float[starFile.RowCount];
            float[] Z = starFile.HasColumn("rlnOriginZAngst") ? starFile.GetColumn("rlnOriginZAngst").Select(v => float.Parse(v, CultureInfo.InvariantCulture)).ToArray() : new float[starFile.RowCount];

            return Helper.Zip(X, Y, Z);
        }

        static void Main(string[] args)
        {

            /*
             * Accepts a single .star and one outdir as input and writes coordinateFiles to outdir for relion picking
             * 
             * */
            if (args.Length != 2)
            {
                Console.WriteLine($"Usage: {System.AppDomain.CurrentDomain.FriendlyName} <starFile with orientations> <outdir for coordinates>");
                return;
            }
            string inFileName = args[0];
            string outDir = args[1];


            Star particleOptics = new Star(inFileName, "optics");

            float pixelSize = float.Parse(particleOptics.GetRowValue(0, "rlnImagePixelSize"), CultureInfo.InvariantCulture);
            Star particleStar = new Star(inFileName, "particles");

            if (!Directory.Exists($@"{outDir}\PredictMonomers\average"))
            {
                Directory.CreateDirectory($@"{outDir}\PredictMonomers\average");
            }
            outDir = $@"{outDir}\PredictMonomers";
            var originalMicrographNames = particleStar.GetRelionMicrographNames().Select(v => "average/" + v.Replace(".tif", ".mrc")).ToArray();

            var imageNames = particleStar.GetColumn("rlnImageName").Select(s =>
            {
                string[] Parts = s.Split('@');
                return Parts[1];
            }).ToArray();
            var uniqueImages = new HashSet<string>(imageNames);
            var angles = particleStar.GetRelionAngles();
            var positions = particleStar.GetRelionCoordinates();
            var origins = GetRelionAngstromOffsets(particleStar);
            var NNMapping = Helper.ArrayOfFunction(idx => -1, positions.Count());
            var distances = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var distanceX = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var distanceY = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var distanceZ = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var matrixDists = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var orientationDistance = Helper.ArrayOfFunction(i => -1.0f, positions.Count());

            Particle[] particles = Helper.ArrayOfFunction(i => new Particle(positions[i] * pixelSize - origins[i], angles[i] * Helper.ToRad, i, originalMicrographNames[i], 0), imageNames.Length);

            List<int>[] particlesInMics = Helper.ArrayOfFunction(i => new List<int>(), uniqueImages.Count());

            for (int idxGlobal = 0, micrographIdx = 0; idxGlobal < imageNames.Length; micrographIdx++)
            {

                var micrographName = imageNames[idxGlobal];

                while (idxGlobal < imageNames.Length && micrographName == imageNames[idxGlobal])
                {
                    particlesInMics[micrographIdx].Add(idxGlobal);
                    particles[idxGlobal].micrographIdx = micrographIdx;
                    idxGlobal++;
                }
            }
            Helper.ForCPU(0, particlesInMics.Count(), 30, null, (micrographIdx, t) =>
            {
                var particlesInMic = particlesInMics[micrographIdx];
                for (int idxInMic = 0; idxInMic < particlesInMic.Count(); idxInMic++)
                {

                    int globalIdx = particlesInMic[idxInMic];
                    Particle particleI = particles[particlesInMic[idxInMic]];
                    float smallestD = float.MaxValue;
                    float smallestO = float.MaxValue;
                    float dX = 0;
                    float dY = 0;
                    float dZ = 0;
                    float matrixDist = 0;

                    int smallestGlobalJdx = -1;
                    for (int jdxInMic = 0; jdxInMic < particlesInMic.Count(); jdxInMic++)
                    {

                        Particle particleJ = particles[particlesInMic[jdxInMic]];
                        if (idxInMic == jdxInMic)
                            continue;
                        float d = (particleI.pos - particleJ.pos).Length();

                        if (d < smallestD)
                        {

                            smallestD = d;
                            smallestGlobalJdx = particleJ.idx;
                            float3 direction = particleJ.pos - particleI.pos;
                            float3 directionRot = particleI.orientation.Transposed() * direction;
                            dX = directionRot.X;
                            dY = directionRot.Y;
                            dZ = directionRot.Z;
                            if(dY > 0)
                            {
                                Matrix3 R = particleJ.orientation * particleI.orientation.Transposed();
                                smallestO = (float)Math.Acos((R.M11 + R.M22 + R.M33 - 1) / 2.0d);  //http://www.boris-belousov.net/2016/12/01/quat-dist/#using-rotation-matrices
                            }
                            else
                            {
                                Matrix3 R = particleI.orientation * particleJ.orientation.Transposed();
                                smallestO = (float)Math.Acos((R.M11 + R.M22 + R.M33 - 1) / 2.0d);  //http://www.boris-belousov.net/2016/12/01/quat-dist/#using-rotation-matrices
                            }
                        }
                    }
                    NNMapping[globalIdx] = smallestGlobalJdx;
                    distances[globalIdx] = smallestD;
                    distanceX[globalIdx] = dX;
                    distanceY[globalIdx] = dY;
                    distanceZ[globalIdx] = dZ;
                    matrixDists[globalIdx] = matrixDist;
                    orientationDistance[globalIdx] = smallestO;
                }
            }, null);

            int newParticleIdx = 1;
            

            float3 TheoDir = new float3(-25, -60, 60);  //very approximate expected orientation
            
            Star filteredParticles = new Star(particleStar.GetColumnNames());
            Helper.ForCPU(0, particlesInMics.Count(), 15, null, (micrographIdx, t) =>
            {
                Star newParticles = new Star(new string[] { "rlnCoordinateX", "rlnCoordinateY" });
                string micrographName = imageNames[particles[particlesInMics[micrographIdx][0]].idx];

                for (int idxInMic = 0; idxInMic < particlesInMics[micrographIdx].Count(); idxInMic++)
                {
                    int globalIdx = particlesInMics[micrographIdx][idxInMic];
                    float3 thisTheoDir = particles[globalIdx].orientation * TheoDir;
                    newParticles.AddRow(new List<string>(new string[] { $"{particles[globalIdx].pos.X / pixelSize}", $"{particles[globalIdx].pos.Y / pixelSize}" }));
                    newParticleIdx++;
                    if (true && !(distances[globalIdx] < 90 && orientationDistance[globalIdx] > 2.9)) //Add new particles based on most probable orientation
                    {
                        float3 newLoc = particles[globalIdx].pos + thisTheoDir;
                        newParticles.AddRow(new List<string>(new string[] { $"{(newLoc.X / pixelSize).ToString(CultureInfo.InvariantCulture)}", $"{(newLoc.Y / pixelSize).ToString(CultureInfo.InvariantCulture)}" }));

                        newParticleIdx++;
                        newLoc = particles[globalIdx].pos - thisTheoDir;
                        newParticles.AddRow(new List<string>(new string[] { $"{(newLoc.X / pixelSize).ToString(CultureInfo.InvariantCulture)}", $"{(newLoc.Y / pixelSize).ToString(CultureInfo.InvariantCulture)}" }));

                        newParticleIdx++;
                    }
                }
                newParticles.Save($@"{outDir}\{micrographName.Replace("particles", "average").Replace("_SARSCoV2_nsp12_net_5.mrcs", "")}_predictedpick.star");

            }, null);
            using (System.IO.StreamWriter file = new System.IO.StreamWriter($@"{outDir}\coords_suffix_predictedpick.star", false))
            {
                file.WriteLine("micrographs.star");
            }
        }

    }
}
