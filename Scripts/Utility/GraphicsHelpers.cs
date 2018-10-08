
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PhysicalSky.Utilities
{
	public class GraphicsHelpers : MonoBehaviour
	{
        public static void DumpTexture(RenderTexture rt, string path)
        {
                Graphics.CopyTexture(rt, i, temp_rt, 0);
                RenderTexture.active = temp_rt;

            }

            temp_rt.Release();

            byte[] bytes = outputTexture.EncodeToEXR(Texture2D.EXRFlags.OutputAsFloat);

            if (!path.ToLower().EndsWith(".exr"))
                path = path + ".exr";

            string directory = System.IO.Path.GetDirectoryName(path);
            if (!System.IO.Directory.Exists(directory))
                System.IO.Directory.CreateDirectory(directory);

            System.IO.File.WriteAllBytes(path, bytes);
        }

		public static void Blit(RenderTexture dest, Material mat, int pass)
		{
			Graphics.Blit(null, dest, mat, pass);
		}

		public static void BlitWithDummy(Material mat, int pass, int width, int height)
		{
			RenderTexture dummy = RenderTexture.GetTemporary(width, height, 0, RenderTextureFormat.R8);
			Blit(dummy, mat, pass);
			RenderTexture.ReleaseTemporary(dummy);
		}

        public static void Blit3D(RenderTargetSetup rtSetup, int depthSlices, Material mat, int pass)
        {
            GL.PushMatrix();
            GL.LoadOrtho();
            for (int i = 0; i < depthSlices; ++i)
            {
                rtSetup.depthSlice = i;
                Graphics.SetRenderTarget(rtSetup);

                float z = (i + 0.5f) / depthSlices;

                mat.SetPass(pass);
                GL.Begin(GL.QUADS);

                GL.TexCoord3(0, 0, z);
                GL.Vertex3(0, 0, 0);
                GL.TexCoord3(1, 0, z);
                GL.Vertex3(1, 0, 0);
                GL.TexCoord3(1, 1, z);
                GL.Vertex3(1, 1, 0);
                GL.TexCoord3(0, 1, z);
                GL.Vertex3(0, 1, 0);

                GL.End();
            }
            GL.PopMatrix();
        }

		public static void Blit3D(RenderTexture dest, Material mat, int pass)
		{
            RenderTargetSetup rtSetup = new RenderTargetSetup(dest.colorBuffer, dest.depthBuffer);

            Blit3D(rtSetup, dest.volumeDepth, mat, pass);
        }

        public static void Blit3D(RenderTexture[] dest, Material mat, int pass)
        {
            RenderBuffer[] rb = new RenderBuffer[dest.Length];
            for (int i = 0; i < dest.Length; i++)
                rb[i] = dest[i].colorBuffer;

            RenderTargetSetup rtSetup = new RenderTargetSetup(rb, dest[0].depthBuffer);

            Blit3D(rtSetup, dest[0].volumeDepth, mat, pass);
        }
    }
}