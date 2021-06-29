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
        static float3[] GetRelionAngstromOffsets(Star starFile)
        {//Compatibility function in case of different relion versions. Always return offsets in Å
            float[] X = starFile.HasColumn("rlnOriginXAngst") ? starFile.GetColumn("rlnOriginXAngst").Select(v => float.Parse(v, CultureInfo.InvariantCulture)).ToArray() : new float[starFile.RowCount];
            float[] Y = starFile.HasColumn("rlnOriginYAngst") ? starFile.GetColumn("rlnOriginYAngst").Select(v => float.Parse(v, CultureInfo.InvariantCulture)).ToArray() : new float[starFile.RowCount];
            float[] Z = starFile.HasColumn("rlnOriginZAngst") ? starFile.GetColumn("rlnOriginZAngst").Select(v => float.Parse(v, CultureInfo.InvariantCulture)).ToArray() : new float[starFile.RowCount];

            return Helper.Zip(X, Y, Z);
        }


        public static float3 MyEulerFromMatrix(Matrix3 a)
        {
            float3 rotated = a * (new float3(0, 0, 1));
            return new float3((float)(Math.PI + Math.Atan2(rotated.Y, rotated.X)), (float)(Math.Acos(rotated.Z / rotated.Length())), 0);
        }

        static void Main(string[] args)
        {
            /*
             * Accepts a single .star as input and writes .txt file containing NN information
             * 
             * */
            if (args.Length != 1)
            {
                Console.WriteLine($"Usage: {System.AppDomain.CurrentDomain.FriendlyName} <starFile with orientations>");
                return;
            }
            string inFileName = args[0];
            var ext = Path.GetExtension(inFileName);
            var prefix = Path.GetFileNameWithoutExtension(inFileName);
            var outdir = Path.GetDirectoryName(inFileName);


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

            //Arrays for saving NN information, in order: idx of partner, distance on micrograph, dx, dy, dz in monomer coordinates, rotation matrix between NN, angle of rotation around eigenvector of relative rotation matrix
            var NNMapping = Helper.ArrayOfFunction(idx => -1, positions.Count());
            var distances = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var distanceX = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var distanceY = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var distanceZ = Helper.ArrayOfFunction(i => -1.0f, positions.Count());
            var relMatrix = Helper.ArrayOfFunction(i => new Matrix3(), positions.Count());
            var orientationDistance = Helper.ArrayOfFunction(i => -1.0f, positions.Count());

            //particle information is saved in the Particle struct for each particle. Thereby, positions are converted to angstrom such that distances will be in Angstrom as well.
            Particle[] particles = Helper.ArrayOfFunction(i => new Particle(positions[i] * pixelSize - origins[i], angles[i] * Helper.ToRad, i, originalMicrographNames[i], 0), imageNames.Length);

            List<int>[] particlesInMics = Helper.ArrayOfFunction(i => new List<int>(), uniqueImages.Count());


            //Set up lookup of micrograph->particles (particlesInMics) and for each particle, set the micrographIdx it originated from
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

            /* find nearest neighbor, parallelize over different micrographs, as they can be handled independently */
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
                    ;
                    Matrix3 relativeMat = new Matrix3();
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
                            dX = directionRot.X;
                            dY = directionRot.Y;
                            dZ = directionRot.Z;

                            //Make an arbitrary choice to make choice of the reference monomer for two NN the same no matter which one is i and j
                            if (dY > 0)
                            {
                                // Calculate relative orientation matrix
                                Matrix3 R = particleJ.orientation * particleI.orientation.Transposed();
                                relativeMat = R;
                                // Relative orientation expressed as a single angle (rotation around the eigenaxis of the rotation matrix)
                                smallestO = (float)Math.Acos((R.M11 + R.M22 + R.M33 - 1) / 2.0d);  //this is the angle of rotation around the eigenvector of R (euler theory of rotation) http://www.boris-belousov.net/2016/12/01/quat-dist/#using-rotation-matrices

                            }
                            else
                            {
                                // Calculate relative orientation matrix
                                Matrix3 R = particleI.orientation * particleJ.orientation.Transposed();
                                relativeMat = R;
                                // Relative orientation expressed as a single angle (rotation around the eigenaxis of the rotation matrix)
                                smallestO = (float)Math.Acos((R.M11 + R.M22 + R.M33 - 1) / 2.0d);  //this is the angle of rotation around the eigenvector of R (euler theory of rotation) http://www.boris-belousov.net/2016/12/01/quat-dist/#using-rotation-matrices
                            }
                        }
                    }
                    NNMapping[globalIdx] = smallestGlobalJdx;
                    distances[globalIdx] = smallestD;
                    distanceX[globalIdx] = dX;
                    distanceY[globalIdx] = dY;
                    distanceZ[globalIdx] = dZ;
                    orientationDistance[globalIdx] = smallestO;
                    relMatrix[globalIdx] = relativeMat;
                }
            }, null);

            //Output relative orientation angles in .star format for compatibility with tools to plot relative orientation
            Star relativeStar = new Star(new String[] { "rlnAngleRot", "rlnAngleTilt", "distance" });
            //This outfile is pandas compatible
            using (System.IO.StreamWriter file = new System.IO.StreamWriter($@"{outdir}\{prefix}.NNOrientations.txt", false))
            {
                file.WriteLine($"particleIdx\tpartnerIdx\tmicrographIdx\tdistance\tdx\tdy\tdz\torientation\tmatrixDist\trot\ttilt\tpsi");

                for (int micrographIdx = 0; micrographIdx < particlesInMics.Count(); micrographIdx++)
                {
                    string micrographName = imageNames[particles[particlesInMics[micrographIdx][0]].idx];

                    for (int idxInMic = 0; idxInMic < particlesInMics[micrographIdx].Count(); idxInMic++)
                    {
                        int globalIdx = particlesInMics[micrographIdx][idxInMic];
                        if (distances[globalIdx] != float.MaxValue)
                        {
                            float3 rotAngles = Matrix3.EulerFromMatrix(relMatrix[globalIdx]);
                            float3 altRotAngles = MyEulerFromMatrix(relMatrix[globalIdx]) * Helper.ToDeg; //Only rot, tilt
                            file.WriteLine($"{globalIdx}\t{NNMapping[globalIdx]}\t{particles[globalIdx].micrographIdx}\t{distances[globalIdx].ToString(CultureInfo.InvariantCulture)}\t{distanceX[globalIdx].ToString(CultureInfo.InvariantCulture)}\t{distanceY[globalIdx].ToString(CultureInfo.InvariantCulture)}\t{distanceZ[globalIdx].ToString(CultureInfo.InvariantCulture)}\t{orientationDistance[globalIdx].ToString(CultureInfo.InvariantCulture)}\t{rotAngles.X.ToString(CultureInfo.InvariantCulture)}\t{rotAngles.Y.ToString(CultureInfo.InvariantCulture)}\t{rotAngles.Z.ToString(CultureInfo.InvariantCulture)}");
                            relativeStar.AddRow(new List<String>(new string[] { $"{altRotAngles.X.ToString(CultureInfo.InvariantCulture)}", $"{altRotAngles.Y.ToString(CultureInfo.InvariantCulture)}", $"{distances[globalIdx].ToString(CultureInfo.InvariantCulture)}" }));

                        }
                        else
                            Console.WriteLine($"WARNIN: No nearest neighbor was found for particle {globalIdx} - probably only one particle on this micrograph");

                    }
                }
            }
            relativeStar.Save($@"{outdir}\{prefix}.NNOrientations.star");

        }

    }
}

