using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PhysicalSky
{
    /// <summary>
    /// Handles the majority of the internal functionality of the AtmosphereModel
    /// </summary>
    public partial class AtmosphereModel
    {
        private class TextureFactory
        {
            // Same as defined in PhysicalSkyCommon
            const int TRANSMITTANCE_TEXTURE_WIDTH = 256;
            const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;
            const int SCATTERING_TEXTURE_R_SIZE = 32;
            const int SCATTERING_TEXTURE_MU_SIZE = 128;
            const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
            const int SCATTERING_TEXTURE_NU_SIZE = 8;
            const int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE;
            static readonly int SCATTERING_TEXTURE_HEIGHT = Supports3DTextures() ? SCATTERING_TEXTURE_MU_SIZE : SCATTERING_TEXTURE_R_SIZE * SCATTERING_TEXTURE_MU_SIZE;
            static readonly int SCATTERING_TEXTURE_DEPTH = Supports3DTextures() ? SCATTERING_TEXTURE_R_SIZE : 1;
            const int IRRADIANCE_TEXTURE_WIDTH = 64;
            const int IRRADIANCE_TEXTURE_HEIGHT = 16;

            public enum Preset
            {
                Transmittance,
                Irradiance,
                Scattering
            }

            public static bool Supports3DTextures()
            {
                return true;
            }

            private static FilterMode GetFilterMode(Preset preset)
            {
                if (preset == Preset.Scattering)
                    return FilterMode.Trilinear;
                else
                    return FilterMode.Bilinear;
            }

            private static RenderTextureDescriptor GetTextureDescriptor(Preset preset)
            {
                RenderTextureDescriptor desc;
                switch (preset)
                {
                    case Preset.Transmittance:
                        desc = new RenderTextureDescriptor(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
                        break;
                    case Preset.Irradiance:
                        desc = new RenderTextureDescriptor(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
                        break;
                    case Preset.Scattering:
                        desc = new RenderTextureDescriptor(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
                        if (Supports3DTextures())
                        {
                            desc.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
                            desc.volumeDepth = SCATTERING_TEXTURE_DEPTH;
                        }
                        break;
                    default:
                        throw new System.NotImplementedException();
                }

                desc.useMipMap = false;
                desc.sRGB = false;
                desc.colorFormat = RenderTextureFormat.ARGBFloat; // TODO: Try changing to half or float
                desc.depthBufferBits = 0;
                return desc;
            }

            public static void CreateRenderTexture(ref RenderTexture tex, Preset preset)
            {
                if (tex == null)
                    tex = new RenderTexture(GetTextureDescriptor(preset));

                if (!tex.IsCreated())
                {
                    tex.wrapMode = TextureWrapMode.Clamp;
                    tex.filterMode = GetFilterMode(preset);
                    tex.Create();
                }
            }

            public static void ReleaseRenderTexture(ref RenderTexture tex)
            {
                if (tex != null)
                {
                    tex.Release();
                }
                tex = null;
            }

            public static RenderTexture CreateTempRenderTexture(Preset preset)
            {
                RenderTexture tex = RenderTexture.GetTemporary(GetTextureDescriptor(preset));
                tex.wrapMode = TextureWrapMode.Clamp;
                tex.filterMode = GetFilterMode(preset);
                return tex;
            }

            public static void ReleaseTempRenderTexture(ref RenderTexture tex)
            {
                RenderTexture.ReleaseTemporary(tex);
                tex = null;
            }
        }
    }
}