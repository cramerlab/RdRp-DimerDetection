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



    class Program
    {
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
            }
            string inFileName = args[0];
            string outDir = args[1];

            if (!Directory.Exists($@"{outDir}\DimerPick\average"))
            {
                Directory.CreateDirectory($@"{outDir}\DimerPick\average");
            }
            outDir = $@"{outDir}\DimerPick";

            Star particleOptics = new Star(inFileName, "optics");

            float pixelSize = float.Parse(particleOptics.GetRowValue(0, "rlnImagePixelSize"), CultureInfo.InvariantCulture);
            Star particleStar = new Star(inFileName, "particles");

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
            var NNMapping = Helper.ArrayOfFunction(idx => new int[] { -1, -1 }, positions.Count());
            var distances = Helper.ArrayOfFunction(i => new float[] { -1f, -1f }, positions.Count());
            var distanceX = Helper.ArrayOfFunction(i => new float[] { -1f, -1f }, positions.Count());
            var distanceY = Helper.ArrayOfFunction(i => new float[] { -1f, -1f }, positions.Count());
            var distanceZ = Helper.ArrayOfFunction(i => new float[] { -1f, -1f }, positions.Count());
            var orientationDistance = Helper.ArrayOfFunction(i => new float[] { -1f, -1f }, positions.Count());

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
                    float[] smallestD = new float[] { float.MaxValue, float.MaxValue };

                    float[] smallestO = new float[] { float.MaxValue, float.MaxValue };
                    float[] dX = new float[] { -1f, -1f };
                    float[] dY = new float[] { -1f, -1f };
                    float[] dZ = new float[] { -1f, -1f };

                    int[] smallestGlobalJdx = new int[] { -1, -1 };
                    for (int jdxInMic = 0; jdxInMic < particlesInMic.Count(); jdxInMic++)
                    {

                        Particle particleJ = particles[particlesInMic[jdxInMic]];
                        if (idxInMic == jdxInMic)
                            continue;
                        float d = (particleI.pos - particleJ.pos).Length();

                        // We save not only the NN (smallestD[0], but also second NN, which is unused in this code)
                        if (d < smallestD[0])
                        {

                            smallestD[0] = d;
                            smallestGlobalJdx[0] = particleJ.idx;
                            Matrix3 R = particleJ.orientation * particleI.orientation.Transposed();
                            smallestO[0] = (float)Math.Acos((R.M11 + R.M22 + R.M33 - 1) / 2.0d);  //http://www.boris-belousov.net/2016/12/01/quat-dist/#using-rotation-matrices
                            float3 direction = particleJ.pos - particleI.pos;
                            /*
                              We make an approximation that the tilt is a perfect 30 degree rotation around the x axis. MEasured distance in y needs to be adjusted
                            */
                            direction.Y = (float)(direction.Y / Math.Cos(30.0 / 360.0 * 2 * Math.PI));

                            float3 directionRot = particleI.orientation.Transposed() * direction;
                            dX[0] = directionRot.X;
                            dY[0] = directionRot.Y;
                            dZ[0] = directionRot.Z;
                        }
                        else if (d < smallestD[1])
                        {
                            smallestD[1] = d;
                            smallestGlobalJdx[1] = particleJ.idx;
                            Matrix3 R = particleJ.orientation * particleI.orientation.Transposed();
                            smallestO[1] = (float)Math.Acos((R.M11 + R.M22 + R.M33 - 1) / 2.0d);  //http://www.boris-belousov.net/2016/12/01/quat-dist/#using-rotation-matrices
                            float3 direction = particleJ.pos - particleI.pos;
                            /*
                              We make an approximation that the tilt is a perfect 30 degree rotation around the x axis. MEasured distance in y needs to be adjusted
                            */
                            direction.Y = (float)(direction.Y / Math.Cos(30.0 / 360.0 * 2 * Math.PI));

                            /*
                              directionRot is the vector that connects one monomer to its NN, expressed in the coordinate system of the monomer.
                              To obtain it, we need to apply the inverse of the monomer pose to the vector that connects both monomers in the
                              micrograph coordinate system
                            */
                            float3 directionRot = particleI.orientation.Transposed() * direction;
                            dX[1] = directionRot.X;
                            dY[1] = directionRot.Y;
                            dZ[1] = directionRot.Z;
                        }
                    }
                    // Save NN information for current (Idx) particle
                    NNMapping[globalIdx] = smallestGlobalJdx;
                    distances[globalIdx] = smallestD;
                    distanceX[globalIdx] = dX;
                    distanceY[globalIdx] = dY;
                    distanceZ[globalIdx] = dZ;
                    orientationDistance[globalIdx] = smallestO;


                }
            }, null);

            // Now we write out all NN that fullfill our dimer criterion
            Helper.ForCPU(0, particlesInMics.Count(), 15, null, (micrographIdx, t) =>
            {
                string micrographName = imageNames[particles[particlesInMics[micrographIdx][0]].idx];
                Star dimerPositions = new Star(new string[] { "rlnCoordinateX", "rlnCoordinateY" });

                for (int idxInMicrograph = 0; idxInMicrograph < particlesInMics[micrographIdx].Count(); idxInMicrograph++)
                {
                    int globalIdx = particlesInMics[micrographIdx][idxInMicrograph];

                    //check if dimer criterion is fullfilled
                    if ((orientationDistance[globalIdx][0] > 2.9) && globalIdx < NNMapping[globalIdx][0] && distances[globalIdx][0] > 20 && distances[globalIdx][0] < 90)
                    {
                        float3 newPos = (particles[globalIdx].pos + particles[NNMapping[globalIdx][0]].pos) / 2;
                        dimerPositions.AddRow(new List<string>(new string[] { $"{(newPos.X / pixelSize).ToString(CultureInfo.InvariantCulture)}", $"{(newPos.Y / pixelSize).ToString(CultureInfo.InvariantCulture)}" }));
                    }

                }

                //Write out non empty .star files
                if (dimerPositions.RowCount > 0)
                    dimerPositions.Save($@"{outDir}\average\{Path.GetFileName(micrographName).Replace(".mrcs", "_dimerPick.star")}");

            }, null);
            using (System.IO.StreamWriter file = new System.IO.StreamWriter($@"{outDir}\coords_suffix_dimerPick.star", false))
            {
                file.WriteLine("micrographs.star");
            }
        }
    }



}


