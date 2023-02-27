using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Vintagestory.API.Common;
using Vintagestory.API.Config;
using Vintagestory.API.Datastructures;
using Vintagestory.API.MathTools;
using Vintagestory.API.Server;
using Vintagestory.API.Util;
using Vintagestory.ServerMods.NoObf;

namespace Vintagestory.ServerMods
{
    public class GenTerra : ModStdWorldGen
    {
        ICoreServerAPI api;
        int regionMapSize;

        private NormalizedSimplexNoise TerrainNoise;

        LandformsWorldProperty landforms;
        Dictionary<int, LerpedWeightedIndex2DMap> LandformMapByRegion = new Dictionary<int, LerpedWeightedIndex2DMap>(10);

        SimplexNoise distort2dx;
        SimplexNoise distort2dz;
        NormalizedSimplexNoise geoUpheavalNoise;
        NormalizedSimplexNoise geoOceanNoise;
        float oceanicityStrengthInv;

        // We generate the whole terrain here so we instantly know the heightmap
        const int lerpHor = TerraGenConfig.lerpHorizontal;
        const int lerpVer = TerraGenConfig.lerpVertical;
        int noiseWidth;
        int paddedNoiseWidth;
        int paddedNoiseHeight;
        int noiseHeight;
        const float lerpDeltaHor = 1f / lerpHor;
        const float lerpDeltaVert = 1f / lerpVer;

        [Obsolete] int previouslyPaddedNoiseWidth;
        [Obsolete] float previouslyWeirdMultiplier;

        double[] noiseTemp;
        float noiseScale;

        float[] terrainThresholdsX0;
        float[] terrainThresholdsX1;
        float[] terrainThresholdsX2;
        float[] terrainThresholdsX3;

        double continentalNoiseOffsetX;
        double continentalNoiseOffsetZ;


        public override bool ShouldLoad(EnumAppSide side)
        {
            return side == EnumAppSide.Server;
        }


        public override double ExecuteOrder()
        {
            return 0;
        }

        public override void StartServerSide(ICoreServerAPI api)
        {
            this.api = api;

            api.Event.ServerRunPhase(EnumServerRunPhase.ModsAndConfigReady, loadGamePre);
            api.Event.InitWorldGenerator(initWorldGen, "standard");
            api.Event.ChunkColumnGeneration(OnChunkColumnGen, EnumWorldGenPass.Terrain, "standard");

            api.Event.SaveGameLoaded += Event_SaveGameLoaded;
        }

        private void Event_SaveGameLoaded()
        {
            continentalNoiseOffsetX = api.WorldManager.SaveGame.GetData<double>("continentalNoiseOffsetX");
            continentalNoiseOffsetZ = api.WorldManager.SaveGame.GetData<double>("continentalNoiseOffsetZ");

            ITreeAttribute worldConfig = api.WorldManager.SaveGame.WorldConfiguration;
            float str = worldConfig.GetString("oceanessStrength", "0").ToFloat(0);

            // oceanicityStrengthInv = 0       => everywhere ocean
            // oceanicityStrengthInv = 0.65    => common oceans
            // oceanicityStrengthInv = 0.85    => rare
            // oceanicityStrengthInv = 1       => no oceans
            oceanicityStrengthInv = GameMath.Clamp(1 - str, 0, 1);
        }

        private void loadGamePre()
        {
            if (api.WorldManager.SaveGame.WorldType != "standard") return;
            
            TerraGenConfig.seaLevel = (int)(0.4313725490196078 * api.WorldManager.MapSizeY);
            api.WorldManager.SetSeaLevel(TerraGenConfig.seaLevel);
        }

        int terrainGenOctaves = 9;
        double[] lerpedAmps;
        double[] lerpedTh;
        float[] distY;

        public void initWorldGen()
        {
            LoadGlobalConfig(api);
            LandformMapByRegion.Clear();

            chunksize = api.WorldManager.ChunkSize;

            regionMapSize = (int)Math.Ceiling((double)api.WorldManager.MapSizeX / api.WorldManager.RegionSize);
            noiseScale = Math.Max(1, api.WorldManager.MapSizeY / 256f);

            terrainGenOctaves = TerraGenConfig.GetTerrainOctaveCount(api.WorldManager.MapSizeY);

            TerrainNoise = NormalizedSimplexNoise.FromDefaultOctaves(
                terrainGenOctaves, 0.0005 / noiseScale, 0.9, api.WorldManager.Seed
            );
            lerpedAmps = new double[terrainGenOctaves];
            lerpedTh = new double[terrainGenOctaves];


            noiseWidth = chunksize / lerpHor;
            noiseHeight = api.WorldManager.MapSizeY / lerpVer;

            paddedNoiseWidth = noiseWidth + 1;
            paddedNoiseHeight = noiseHeight + 1;

            // Until v1.18.0, a few key calculations used paddedNoiseWidth when they should have used noiseWidth.
            if (GameVersion.IsAtLeastVersion(api.WorldManager.SaveGame.CreatedGameVersion, "1.18.0-pre.1"))
            {
                previouslyPaddedNoiseWidth = noiseWidth;
                previouslyWeirdMultiplier = 1.0f;
            }
            else
            {
                previouslyPaddedNoiseWidth = paddedNoiseWidth;
                previouslyWeirdMultiplier = (float)lerpHor / chunksize + 1.0f;
            }

            noiseTemp = new double[paddedNoiseWidth * paddedNoiseWidth * paddedNoiseHeight];
            distY = new float[paddedNoiseWidth * paddedNoiseWidth];

            distort2dx = new SimplexNoise(
                new double[] { 55, 40, 30, 10 }, 
                scaleAdjustedFreqs(new double[] { 1 / 5.0, 1 / 2.50, 1 / 1.250, 1 / 0.65 }, noiseScale), 
                api.World.Seed + 9876 + 0
            );
            distort2dz = new SimplexNoise(
                new double[] { 55, 40, 30, 10 },
                scaleAdjustedFreqs(new double[] { 1 / 5.0, 1 / 2.50, 1 / 1.250, 1 / 0.65 }, noiseScale),
                api.World.Seed + 9876 + 2
            );
            geoUpheavalNoise = new NormalizedSimplexNoise(
                new double[] { 55, 40, 30, 10, 2 },
                scaleAdjustedFreqs(new double[] { 1 / 5.0 / 1.1, 1 / 2.50 / 1.1, 1 / 1.250 / 1.1, 1 / 0.65 / 1.1, 4 }, noiseScale), 
                api.World.Seed + 9876 + 1
            );

            geoOceanNoise = NormalizedSimplexNoise.FromDefaultOctaves(6, 1 / 60.0, 0.8, api.World.Seed + 9856 + 1);

            // Find a non-oceanic spot
            if (api.WorldManager.SaveGame.IsNew && oceanicityStrengthInv < 1)
            {
                var rnd = new Random(api.World.Seed + 2837986);
                BlockPos pos = new BlockPos(api.WorldManager.MapSizeX / 2, 0, api.WorldManager.MapSizeZ / 2);
                int tries = 0;

                while (tries++ < 4000)
                {
                    continentalNoiseOffsetX = (10 * tries) / 400.0 * (1 - 2 * rnd.Next(2));
                    continentalNoiseOffsetZ = (10 * tries) / 400.0 * (1 - 2 * rnd.Next(2));

                    var noiseVal = (int)Math.Max(0, 255 * Math.Min(1, 2 * geoOceanNoise.Noise(
                        continentalNoiseOffsetX + pos.X / 400.0,
                        continentalNoiseOffsetZ + pos.Z / 400.0
                    ) - oceanicityStrengthInv));
                    
                    if (noiseVal <= 0)
                    {
                        api.WorldManager.SaveGame.StoreData<double>("continentalNoiseOffsetX", continentalNoiseOffsetX);
                        api.WorldManager.SaveGame.StoreData<double>("continentalNoiseOffsetZ", continentalNoiseOffsetZ);
                        break;
                    }
                }
            }

            terrainThresholdsX0 = new float[api.WorldManager.MapSizeY];
            terrainThresholdsX1 = new float[api.WorldManager.MapSizeY];
            terrainThresholdsX2 = new float[api.WorldManager.MapSizeY];
            terrainThresholdsX3 = new float[api.WorldManager.MapSizeY];

            api.Logger.VerboseDebug("Initialised GenTerra");
        }

        private double[] scaleAdjustedFreqs(double[] vs, float horizontalScale)
        {
            for (int i = 0; i < vs.Length; i++)
            {
                vs[i] /= horizontalScale;
            }

            return vs;
        }


        

        private void OnChunkColumnGen(IServerChunk[] chunks, int chunkX, int chunkZ, ITreeAttribute chunkGenParams = null)
        {
            landforms = NoiseLandforms.landforms;
            IMapChunk mapchunk = chunks[0].MapChunk;
            int chunksize = this.chunksize;

            int climateUpLeft;
            int climateUpRight;
            int climateBotLeft;
            int climateBotRight;

            int upheavalMapUpLeft = 0;
            int upheavalMapUpRight = 0;
            int upheavalMapBotLeft = 0;
            int upheavalMapBotRight = 0;

            IntDataMap2D climateMap = chunks[0].MapChunk.MapRegion.ClimateMap;
            int regionChunkSize = api.WorldManager.RegionSize / chunksize;
            float cfac = (float)climateMap.InnerSize / regionChunkSize;
            int rlX = chunkX % regionChunkSize;
            int rlZ = chunkZ % regionChunkSize;

            climateUpLeft = climateMap.GetUnpaddedInt((int)(rlX * cfac), (int)(rlZ * cfac));
            climateUpRight = climateMap.GetUnpaddedInt((int)(rlX * cfac + cfac), (int)(rlZ * cfac));
            climateBotLeft = climateMap.GetUnpaddedInt((int)(rlX * cfac), (int)(rlZ * cfac + cfac));
            climateBotRight = climateMap.GetUnpaddedInt((int)(rlX * cfac + cfac), (int)(rlZ * cfac + cfac));

            IntDataMap2D upheavalMap = chunks[0].MapChunk.MapRegion.UpheavelMap;
            if (upheavalMap != null)
            {
                float ufac = (float)upheavalMap.InnerSize / regionChunkSize;
                upheavalMapUpLeft = upheavalMap.GetUnpaddedInt((int)(rlX * ufac), (int)(rlZ * ufac));
                upheavalMapUpRight = upheavalMap.GetUnpaddedInt((int)(rlX * ufac + ufac), (int)(rlZ * ufac));
                upheavalMapBotLeft = upheavalMap.GetUnpaddedInt((int)(rlX * ufac), (int)(rlZ * ufac + ufac));
                upheavalMapBotRight = upheavalMap.GetUnpaddedInt((int)(rlX * ufac + ufac), (int)(rlZ * ufac + ufac));
            }

            int waterID = GlobalConfig.waterBlockId;
            int rockID = GlobalConfig.defaultRockId;


            IntDataMap2D landformMap = mapchunk.MapRegion.LandformMap;
            // Amount of pixels for each chunk (probably 1, 2, or 4) in the land form map
            float chunkPixelSize = landformMap.InnerSize / regionChunkSize;
            // Adjusted lerp for the noiseWidth
            float chunkPixelStep = chunkPixelSize / noiseWidth;
            // Start coordinates for the chunk in the region map
            float baseX = (chunkX % regionChunkSize) * chunkPixelSize;
            float baseZ = (chunkZ % regionChunkSize) * chunkPixelSize;
            

            LerpedWeightedIndex2DMap landLerpMap = GetOrLoadLerpedLandformMap(chunks[0].MapChunk, chunkX / regionChunkSize, chunkZ / regionChunkSize);

            // Terrain octaves
            double[] octNoiseX0, octNoiseX1, octNoiseX2, octNoiseX3;
            double[] octThX0, octThX1, octThX2, octThX3;

            // Correct chunk boundary errors prior to 1.18.0
            chunkPixelSize *= previouslyWeirdMultiplier;

            GetInterpolatedOctaves(landLerpMap[baseX, baseZ], out octNoiseX0, out octThX0);
            GetInterpolatedOctaves(landLerpMap[baseX + chunkPixelSize, baseZ], out octNoiseX1, out octThX1);
            GetInterpolatedOctaves(landLerpMap[baseX, baseZ + chunkPixelSize], out octNoiseX2, out octThX2);
            GetInterpolatedOctaves(landLerpMap[baseX + chunkPixelSize, baseZ + chunkPixelSize], out octNoiseX3, out octThX3);


            double[] terrainNoise3d = GetTerrainNoise3D(octNoiseX0, octNoiseX1, octNoiseX2, octNoiseX3, octThX0, octThX1, octThX2, octThX3, chunkX * noiseWidth, 0, chunkZ * noiseWidth);

            // Store heightmap in the map chunk
            ushort[] rainheightmap = chunks[0].MapChunk.RainHeightMap;
            ushort[] terrainheightmap = chunks[0].MapChunk.WorldGenTerrainHeightMap;
            

            // Terrain thresholds
            double tnoiseY0;
            double tnoiseY1;
            double tnoiseY2;
            double tnoiseY3;
            double tnoiseGainY0;
            double tnoiseGainY1;
            double tnoiseGainY2;
            double tnoiseGainY3;


            double thNoiseX0;
            double thNoiseX1;
            double thNoiseGainX0;
            double thNoiseGainX1;
            double thNoiseGainZ0;
            double thNoiseZ0;

            

            

            int mapsizeY = api.WorldManager.MapSizeY;
            int mapsizeYm1 = api.WorldManager.MapSizeY - 1;
            int taperThreshold = (int)(mapsizeY * 0.9f);
            double geoUpheavalAmplitude = 255;

            for (int xN = 0; xN < paddedNoiseWidth; xN++)
            {
                for (int zN = 0; zN < paddedNoiseWidth; zN++)
                {
                    VectorXZ distGeo = DistortionNoise(chunkX * chunksize + xN * lerpHor, chunkZ * chunksize + zN * lerpHor);
                    distGeo = ApplyDistortionThreshold(distGeo * 10.0, 10.0);

                    float upheavalStrength = GameMath.BiLerp(upheavalMapUpLeft, upheavalMapUpRight, upheavalMapBotLeft, upheavalMapBotRight,
                        xN * (1.0f / noiseWidth), zN * (1.0f / noiseWidth));
                    float distYHere = ComputeOceanAndUpheavalDistY(upheavalStrength,
                                chunkX * chunksize + xN * lerpHor, chunkZ * chunksize + zN * lerpHor,
                                distGeo);

                    distY[NoiseIndex2d(xN, zN)] = distYHere;
                }
            }

            for (int xN = 0; xN < noiseWidth; xN++)
            {
                for (int zN = 0; zN < noiseWidth; zN++)
                {
                    // Landform terrain thresholds
                    LoadInterpolatedThresholds(landLerpMap[baseX + xN * chunkPixelStep, baseZ + zN * chunkPixelStep], terrainThresholdsX0);
                    LoadInterpolatedThresholds(landLerpMap[baseX + (xN+1) * chunkPixelStep, baseZ + zN * chunkPixelStep], terrainThresholdsX1);
                    LoadInterpolatedThresholds(landLerpMap[baseX + xN * chunkPixelStep, baseZ + (zN+1) * chunkPixelStep], terrainThresholdsX2);
                    LoadInterpolatedThresholds(landLerpMap[baseX + (xN+1) * chunkPixelStep, baseZ + (zN+1) * chunkPixelStep], terrainThresholdsX3);

                    // Oceanicity/Upheaval Y displacements
                    float distY0 = distY[NoiseIndex2d(xN, zN)];
                    float distY1 = distY[NoiseIndex2d(xN + 1, zN)];
                    float distY2 = distY[NoiseIndex2d(xN, zN + 1)];
                    float distY3 = distY[NoiseIndex2d(xN + 1, zN + 1)];

                    for (int yN = 0; yN < noiseHeight; yN++)
                    {
                        // Terrain noise
                        tnoiseY0 = terrainNoise3d[NoiseIndex3d(xN, yN, zN)];
                        tnoiseY1 = terrainNoise3d[NoiseIndex3d(xN, yN, zN + 1)];
                        tnoiseY2 = terrainNoise3d[NoiseIndex3d(xN + 1, yN, zN)];
                        tnoiseY3 = terrainNoise3d[NoiseIndex3d(xN + 1, yN, zN + 1)];

                        tnoiseGainY0 = (terrainNoise3d[NoiseIndex3d(xN, yN + 1, zN)] - tnoiseY0) * lerpDeltaVert;
                        tnoiseGainY1 = (terrainNoise3d[NoiseIndex3d(xN, yN + 1, zN + 1)] - tnoiseY1) * lerpDeltaVert;
                        tnoiseGainY2 = (terrainNoise3d[NoiseIndex3d(xN + 1, yN + 1, zN)] - tnoiseY2) * lerpDeltaVert;
                        tnoiseGainY3 = (terrainNoise3d[NoiseIndex3d(xN + 1, yN + 1, zN + 1)] - tnoiseY3) * lerpDeltaVert;

                        for (int y = 0; y < lerpVer; y++)
                        {
                            int posY = yN * lerpVer + y;
                            int chunkY = posY / chunksize;
                            int localY = posY % chunksize;
                            IChunkBlocks chunkBlockData = chunks[chunkY].Data;

                            if (posY == 0)
                            {
                                int chunkIndex = ChunkIndex3d(xN * lerpHor, localY, zN * lerpHor);
                                chunkBlockData.SetBlockBulk(chunkIndex, lerpHor, lerpHor, GlobalConfig.mantleBlockId);
                            }
                            else
                            {
                                // For Terrain noise 
                                double tnoiseX0 = tnoiseY0;
                                double tnoiseX1 = tnoiseY1;

                                double tnoiseGainX0 = (tnoiseY2 - tnoiseY0) * lerpDeltaHor;
                                double tnoiseGainX1 = (tnoiseY3 - tnoiseY1) * lerpDeltaHor;

                                thNoiseX0 = SampleDisplacedYThreshold(posY + distY0, mapsizeY - 2, terrainThresholdsX0);
                                thNoiseX1 = SampleDisplacedYThreshold(posY + distY2, mapsizeY - 2, terrainThresholdsX2);
                                double thNoiseY2 = SampleDisplacedYThreshold(posY + distY1, mapsizeY - 2, terrainThresholdsX1);
                                double thNoiseY3 = SampleDisplacedYThreshold(posY + distY3, mapsizeY - 2, terrainThresholdsX3);

                                if (posY >= TerraGenConfig.seaLevel && Math.Max(Math.Max(tnoiseY0, tnoiseY1), Math.Max(tnoiseY2, tnoiseY3)) <= Math.Min(Math.Min(thNoiseX0, thNoiseX1), Math.Min(thNoiseY2, thNoiseY3)))
                                {
                                    // Nothing to do: whole slice is air
                                }
                                else
                                {
                                    thNoiseGainX0 = (thNoiseY2 - thNoiseX0) * lerpDeltaHor;
                                    thNoiseGainX1 = (thNoiseY3 - thNoiseX1) * lerpDeltaHor;

                                    for (int x = 0; x < lerpHor; x++)
                                    {
                                        // For terrain noise
                                        double tnoiseZ0 = tnoiseX0;
                                        double tnoiseGainZ0 = (tnoiseX1 - tnoiseX0) * lerpDeltaHor;

                                        // Landform
                                        thNoiseZ0 = thNoiseX0;
                                        thNoiseGainZ0 = (thNoiseX1 - thNoiseX0) * lerpDeltaHor;

                                        int lX = xN * lerpHor + x;

                                        for (int z = 0; z < lerpHor; z++)
                                        {                                          
                                            int lZ = zN * lerpHor + z;

                                            // Geo Upheaval modifier for threshold
                                            double distYHere = GameMath.BiLerp(distY0, distY2, distY1, distY3, x * lerpDeltaHor, z * lerpDeltaHor);
                                            double geoUpheavalTaper = ComputeGeoUpheavalTaper(posY, distYHere, taperThreshold, geoUpheavalAmplitude, mapsizeY);

                                            if (tnoiseZ0 > thNoiseZ0 + geoUpheavalTaper)
                                            {
                                                int mapIndex = ChunkIndex2d(lX, lZ);
                                                terrainheightmap[mapIndex] = (ushort)posY;
                                                rainheightmap[mapIndex] = (ushort)posY;

                                                chunkBlockData[ChunkIndex3d(lX, localY, lZ)] = rockID;
                                            }
                                            else if (posY < TerraGenConfig.seaLevel)
                                            {
                                                int mapIndex = ChunkIndex2d(lX, lZ);
                                                rainheightmap[mapIndex] = (ushort)posY;

                                                int blockId;
                                                if (posY == TerraGenConfig.seaLevel - 1)
                                                {
                                                    int temp = (GameMath.BiLerpRgbColor(((float)lX) / chunksize, ((float)lZ) / chunksize, climateUpLeft, climateUpRight, climateBotLeft, climateBotRight) >> 16) & 0xff;
                                                    float distort = (float)distort2dx.Noise(chunkX * chunksize + lX, chunkZ * chunksize + lZ) / 20f;
                                                    float tempf = TerraGenConfig.GetScaledAdjustedTemperatureFloat(temp, 0) + distort;

                                                    blockId = (tempf < TerraGenConfig.WaterFreezingTempOnGen) ? GlobalConfig.lakeIceBlockId : waterID;
                                                }
                                                else
                                                {
                                                    blockId = waterID;
                                                }

                                                chunkBlockData.SetFluid(ChunkIndex3d(lX, localY, lZ), blockId);
                                            }

                                            tnoiseZ0 += tnoiseGainZ0;
                                            thNoiseZ0 += thNoiseGainZ0;
                                        }

                                        tnoiseX0 += tnoiseGainX0;
                                        thNoiseX0 += thNoiseGainX0;

                                        tnoiseX1 += tnoiseGainX1;
                                        thNoiseX1 += thNoiseGainX1;
                                    }
                                }
                            }

                            tnoiseY0 += tnoiseGainY0;
                            tnoiseY1 += tnoiseGainY1;
                            tnoiseY2 += tnoiseGainY2;
                            tnoiseY3 += tnoiseGainY3;
                        }
                    }
                }
            }

            ushort ymax = 0;
            for (int i = 0; i < rainheightmap.Length; i++)
            {
                ymax = Math.Max(ymax, rainheightmap[i]);
            }
            chunks[0].MapChunk.YMax = ymax;
        }


        
        LerpedWeightedIndex2DMap GetOrLoadLerpedLandformMap(IMapChunk mapchunk, int regionX, int regionZ)
        {
            LerpedWeightedIndex2DMap map;
            // 1. Load?
            LandformMapByRegion.TryGetValue(regionZ * regionMapSize + regionX, out map);
            if (map != null) return map;

            IntDataMap2D lmap = mapchunk.MapRegion.LandformMap;
            // 2. Create
            map = LandformMapByRegion[regionZ * regionMapSize + regionX] 
                = new LerpedWeightedIndex2DMap(lmap.Data, lmap.Size, TerraGenConfig.landFormSmoothingRadius, lmap.TopLeftPadding, lmap.BottomRightPadding);

            return map;
        }


        // Can be called only once per x/z coordinate to get a list of all thresholds for this column
        private void LoadInterpolatedThresholds(WeightedIndex[] indices, float[] values)
        {
            for (int y = 0; y < values.Length; y++)
            {
                float threshold = 0;
                for (int i = 0; i < indices.Length; i++)
                {
                    threshold += landforms.LandFormsByIndex[indices[i].Index].TerrainYThresholds[y] * indices[i].Weight;
                }

                values[y] = threshold;
            }
        }




        private void GetInterpolatedOctaves(WeightedIndex[] indices, out double[] amps, out double[] thresholds)
        {
            amps = new double[terrainGenOctaves];
            thresholds = new double[terrainGenOctaves];

            for (int octave = 0; octave < terrainGenOctaves; octave++)
            {
                double amplitude = 0;
                double threshold = 0;
                for (int i = 0; i < indices.Length; i++)
                {
                    LandformVariant l = landforms.LandFormsByIndex[indices[i].Index];
                    amplitude += l.TerrainOctaves[octave] * indices[i].Weight;
                    threshold += l.TerrainOctaveThresholds[octave] * indices[i].Weight;
                }

                amps[octave] = amplitude;
                thresholds[octave] = threshold;
            }    
        }


        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        float SampleDisplacedYThreshold(float distortedPosY, int mapSizeYm2, float[] thresholds)
        {
            int distortedPosYBase = (int)Math.Floor(distortedPosY);
            int yBase = GameMath.Clamp(distortedPosYBase, 0, mapSizeYm2);
            float ySlide = distortedPosY - distortedPosYBase;
            return GameMath.Lerp(thresholds[yBase], thresholds[yBase + 1], ySlide);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        float ComputeOceanAndUpheavalDistY(float upheavalStrength, double worldX, double worldZ, VectorXZ distGeo)
        {
            const float OCEANICITY_STRENGTH_MODIFIER = 1.0f / 255.0f;
            float oceanicity = oceanicityStrengthInv >= 1 ? -255 : 255 * Math.Min(1, 2 * (float)geoOceanNoise.Noise(
                continentalNoiseOffsetX + worldX / 400.0,
                continentalNoiseOffsetZ + worldZ / 400.0
            ) - oceanicityStrengthInv);
            if (oceanicity < 0)
            {
                float strength = Math.Min(1, -oceanicity * OCEANICITY_STRENGTH_MODIFIER);
                float upheavalNoiseValue = (float)geoUpheavalNoise.Noise((worldX + distGeo.X) / 400.0, (worldZ + distGeo.Z) / 400.0);
                float upheavalMultiplier = Math.Min(0, 0.5f - upheavalNoiseValue);
                return strength * upheavalStrength * upheavalMultiplier;
            }
            else
            {
                return oceanicity;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        double ComputeGeoUpheavalTaper(double posY, double distY, double taperThreshold, double geoUpheavalAmplitude, double mapSizeY)
        {
            const double AMPLITUDE_MODIFIER = 40.0;
            if (posY > taperThreshold && distY < -2)
            {
                double upheavalAmount = GameMath.Clamp(-distY, posY - mapSizeY, posY);
                double ceilingDelta = posY - taperThreshold;
                return (ceilingDelta * upheavalAmount) / (AMPLITUDE_MODIFIER * geoUpheavalAmplitude);
            }
            return 0;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        VectorXZ DistortionNoise(double worldX, double worldZ)
        {
            double noiseX = worldX / 400.0;
            double noiseZ = worldZ / 400.0;
            double distX = distort2dx.Noise(noiseX, noiseZ);
            double distZ = distort2dz.Noise(noiseX, noiseZ);
            return new VectorXZ { X = distX, Z = distZ };
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        VectorXZ ApplyDistortionThreshold(VectorXZ dist, double threshold)
        {
            dist.X = dist.X > 0 ? Math.Max(0, dist.X - threshold) : Math.Min(0, dist.X + threshold);
            dist.Z = dist.Z > 0 ? Math.Max(0, dist.Z - threshold) : Math.Min(0, dist.Z + threshold);
            return dist;
        }


        double[] GetTerrainNoise3D(double[] octX0, double[] octX1, double[] octX2, double[] octX3, double[] octThX0, double[] octThX1, double[] octThX2, double[] octThX3, int xPos, int yPos, int zPos)
        {
            for (int x = 0; x < paddedNoiseWidth; x++)
            {
                for (int z = 0; z < paddedNoiseWidth; z++)
                {
                    for (int i = 0; i < terrainGenOctaves; i++)
                    {
                        lerpedAmps[i] = GameMath.BiLerp(octX0[i], octX1[i], octX2[i], octX3[i], (double)x / previouslyPaddedNoiseWidth, (double)z / previouslyPaddedNoiseWidth);
                        lerpedTh[i] = GameMath.BiLerp(octThX0[i], octThX1[i], octThX2[i], octThX3[i], (double)x / previouslyPaddedNoiseWidth, (double)z / previouslyPaddedNoiseWidth);
                    }

                    int gridX = xPos + x;
                    int gridZ = zPos + z;
                    int worldX = gridX * lerpHor;
                    int worldZ = gridZ * lerpHor;

                    VectorXZ distTerrain = DistortionNoise(worldX, worldZ);
                    distTerrain = ApplyDistortionThreshold(distTerrain * 4.0, 40.0);

                    for (int y = 0; y < paddedNoiseHeight; y++)
                    {
                        noiseTemp[NoiseIndex3d(x, y, z)] = TerrainNoise.Noise(
                            worldX + distTerrain.X,
                            (yPos + y) * (lerpVer * 0.5 / TerraGenConfig.terrainNoiseVerticalScale),
                            worldZ + distTerrain.Z,
                            lerpedAmps,
                            lerpedTh
                        );
                    }
                }
            }

            return noiseTemp;
        }



        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int ChunkIndex3d(int x, int y, int z)
        {
            return (y * chunksize + z) * chunksize + x;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int ChunkIndex2d(int x, int z)
        {
            return z * chunksize + x;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int NoiseIndex2d(int x, int z)
        {
            return z * paddedNoiseWidth + x;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int NoiseIndex3d(int x, int y, int z)
        {
            return (y * paddedNoiseWidth + z) * paddedNoiseWidth + x;
        }


        struct VectorXZ
        {
            public double X, Z;
            public static VectorXZ operator *(VectorXZ a, double b) => new VectorXZ { X = a.X * b, Z = a.Z * b };
        }
    }
}
