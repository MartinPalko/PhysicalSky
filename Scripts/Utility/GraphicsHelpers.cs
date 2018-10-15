using UnityEngine;

namespace PhysicalSky.Utilities
{
	public class GraphicsHelpers : MonoBehaviour
	{
        private static void Blit3D(RenderTargetSetup rtSetup, int depthSlices, Material mat, int pass)
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

        public static void DumpTexture(RenderTexture rt, string path)
        {
            // Warning: doesn't seem to actually work with 3d textures
            RenderTexture temp_rt = new RenderTexture(rt.width, rt.height, 0, rt.format);
            Texture2D outputTexture = new Texture2D(rt.width, rt.height * rt.volumeDepth, TextureFormat.RGBAFloat, false, true);

            for (int i = 0; i < rt.volumeDepth; i++)
            {
                // CopyTexture seems to always use the first slice anyway, instead of i.
                Graphics.CopyTexture(rt, i, temp_rt, 0);
                RenderTexture.active = temp_rt;

                // Doesn't work at all
                //RenderTargetSetup rtSetup = new RenderTargetSetup(rt.colorBuffer, rt.depthBuffer);
                //rtSetup.depthSlice = i;
                //Graphics.SetRenderTarget(rtSetup);

                outputTexture.ReadPixels(new Rect(0, 0, rt.width, rt.height), 0, rt.height * i, false);
            }

            RenderTexture.active = null;
            temp_rt.Release();

            byte[] bytes = outputTexture.EncodeToEXR(Texture2D.EXRFlags.OutputAsFloat);

            if (!path.ToLower().EndsWith(".exr"))
                path = path + ".exr";

            string directory = System.IO.Path.GetDirectoryName(path);
            if (!System.IO.Directory.Exists(directory))
                System.IO.Directory.CreateDirectory(directory);

            System.IO.File.WriteAllBytes(path, bytes);
        }
    }
}