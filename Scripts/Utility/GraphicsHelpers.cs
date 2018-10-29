using System.Collections.Generic;
using UnityEngine;

namespace PhysicalSky.Utilities
{
	public static class GraphicsHelpers
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

        public static TextureFormat TextureFormatFromRenderTextureFormat(RenderTextureFormat rtf)
        {
            switch(rtf)
            {
                case RenderTextureFormat.ARGBHalf:
                    return TextureFormat.RGBAHalf;
                case RenderTextureFormat.ARGBFloat:
                    return TextureFormat.RGBAFloat;
                default:
                    throw new System.NotImplementedException();
            }
        }

        public static Texture2D Texture2DFromRenderTexture(RenderTexture rt)
        {
            Debug.Assert(rt.volumeDepth == 1);

            Texture2D result = new Texture2D(
                    rt.width,
                    rt.height,
                    TextureFormatFromRenderTextureFormat(rt.format),
                    rt.useMipMap);

            Graphics.SetRenderTarget(rt);
            result.ReadPixels(new Rect(0, 0, rt.width, rt.height), 0, 0, false);
            result.wrapMode = rt.wrapMode;
            result.filterMode = rt.filterMode;
            return result;
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

        public interface IMaterialProperties
        {
            void SetColor(int nameID, Color value);
            void SetColor(string name, Color value);
            void SetFloat(int nameID, float value);
            void SetFloat(string name, float value);
            void SetFloatArray(int nameID, float[] values);
            void SetFloatArray(string name, List<float> values);
            void SetFloatArray(string name, float[] values);
            void SetFloatArray(int nameID, List<float> values);
            void SetInt(int nameID, int value);
            void SetInt(string name, int value);
            void SetMatrix(int nameID, Matrix4x4 value);
            void SetMatrix(string name, Matrix4x4 value);
            void SetMatrixArray(string name, Matrix4x4[] values);
            void SetMatrixArray(int nameID, Matrix4x4[] values);
            void SetMatrixArray(int nameID, List<Matrix4x4> values);
            void SetMatrixArray(string name, List<Matrix4x4> values);
            void SetTexture(int nameID, Texture value);
            void SetTexture(string name, Texture value);
            void SetVector(int nameID, Vector4 value);
            void SetVector(string name, Vector4 value);
            void SetVectorArray(string name, Vector4[] values);
            void SetVectorArray(string name, List<Vector4> values);
            void SetVectorArray(int nameID, Vector4[] values);
            void SetVectorArray(int nameID, List<Vector4> values);
        }

        private class MaterialImpl : IMaterialProperties
        {
            private Material m_material;

            public MaterialImpl(Material material)
            {
                m_material = material;
            }

            public void SetColor(int nameID, Color value) { m_material.SetColor(nameID, value); }
            public void SetColor(string name, Color value) { m_material.SetColor(name, value); }
            public void SetFloat(int nameID, float value) { m_material.SetFloat(nameID, value); }
            public void SetFloat(string name, float value) { m_material.SetFloat(name, value); }
            public void SetFloatArray(int nameID, float[] values) { m_material.SetFloatArray(nameID, values); }
            public void SetFloatArray(string name, List<float> values) { m_material.SetFloatArray(name, values); }
            public void SetFloatArray(string name, float[] values) { m_material.SetFloatArray(name, values); }
            public void SetFloatArray(int nameID, List<float> values) { m_material.SetFloatArray(nameID, values); }
            public void SetInt(int nameID, int value) { m_material.SetInt(nameID, value); }
            public void SetInt(string name, int value) { m_material.SetInt(name, value); }
            public void SetMatrix(int nameID, Matrix4x4 value) { m_material.SetMatrix(nameID, value); }
            public void SetMatrix(string name, Matrix4x4 value) { m_material.SetMatrix(name, value); }
            public void SetMatrixArray(string name, Matrix4x4[] values) { m_material.SetMatrixArray(name, values); }
            public void SetMatrixArray(int nameID, Matrix4x4[] values) { m_material.SetMatrixArray(nameID, values); }
            public void SetMatrixArray(int nameID, List<Matrix4x4> values) { m_material.SetMatrixArray(nameID, values); }
            public void SetMatrixArray(string name, List<Matrix4x4> values) { m_material.SetMatrixArray(name, values); }
            public void SetTexture(int nameID, Texture value) { m_material.SetTexture(nameID, value); }
            public void SetTexture(string name, Texture value) { m_material.SetTexture(name, value); }
            public void SetVector(int nameID, Vector4 value) { m_material.SetVector(nameID, value); }
            public void SetVector(string name, Vector4 value) { m_material.SetVector(name, value); }
            public void SetVectorArray(string name, Vector4[] values) { m_material.SetVectorArray(name, values); }
            public void SetVectorArray(string name, List<Vector4> values) { m_material.SetVectorArray(name, values); }
            public void SetVectorArray(int nameID, Vector4[] values) { m_material.SetVectorArray(nameID, values); }
            public void SetVectorArray(int nameID, List<Vector4> values) { m_material.SetVectorArray(nameID, values); }
        }

        public static IMaterialProperties ToMaterialPropertyInterface(this Material m)
        {
            return new MaterialImpl(m);
        }

        private class MaterialPropertyBlockImpl : IMaterialProperties
        {
            private MaterialPropertyBlock m_block;

            public MaterialPropertyBlockImpl(MaterialPropertyBlock block)
            {
                m_block = block;
            }

            public void SetColor(int nameID, Color value) { m_block.SetColor(nameID, value); }
            public void SetColor(string name, Color value) { m_block.SetColor(name, value); }
            public void SetFloat(int nameID, float value) { m_block.SetFloat(nameID, value); }
            public void SetFloat(string name, float value) { m_block.SetFloat(name, value); }
            public void SetFloatArray(int nameID, float[] values) { m_block.SetFloatArray(nameID, values); }
            public void SetFloatArray(string name, List<float> values) { m_block.SetFloatArray(name, values); }
            public void SetFloatArray(string name, float[] values) { m_block.SetFloatArray(name, values); }
            public void SetFloatArray(int nameID, List<float> values) { m_block.SetFloatArray(nameID, values); }
            public void SetInt(int nameID, int value) { m_block.SetInt(nameID, value); }
            public void SetInt(string name, int value) { m_block.SetInt(name, value); }
            public void SetMatrix(int nameID, Matrix4x4 value) { m_block.SetMatrix(nameID, value); }
            public void SetMatrix(string name, Matrix4x4 value) { m_block.SetMatrix(name, value); }
            public void SetMatrixArray(string name, Matrix4x4[] values) { m_block.SetMatrixArray(name, values); }
            public void SetMatrixArray(int nameID, Matrix4x4[] values) { m_block.SetMatrixArray(nameID, values); }
            public void SetMatrixArray(int nameID, List<Matrix4x4> values) { m_block.SetMatrixArray(nameID, values); }
            public void SetMatrixArray(string name, List<Matrix4x4> values) { m_block.SetMatrixArray(name, values); }
            public void SetTexture(int nameID, Texture value) { m_block.SetTexture(nameID, value); }
            public void SetTexture(string name, Texture value) { m_block.SetTexture(name, value); }
            public void SetVector(int nameID, Vector4 value) { m_block.SetVector(nameID, value); }
            public void SetVector(string name, Vector4 value) { m_block.SetVector(name, value); }
            public void SetVectorArray(string name, Vector4[] values) { m_block.SetVectorArray(name, values); }
            public void SetVectorArray(string name, List<Vector4> values) { m_block.SetVectorArray(name, values); }
            public void SetVectorArray(int nameID, Vector4[] values) { m_block.SetVectorArray(nameID, values); }
            public void SetVectorArray(int nameID, List<Vector4> values) { m_block.SetVectorArray(nameID, values); }
        }

        public static IMaterialProperties ToMaterialPropertyInterface(this MaterialPropertyBlock m)
        {
            return new MaterialPropertyBlockImpl(m);
        }

        static void SetProperties(IMaterialProperties m)
        {
            m.SetVector("Dummy", Vector4.zero);
        }

        static void Test(Material m)
        {
            SetProperties(m.ToMaterialPropertyInterface());
        }

            class UserSurrogate : User
    {
        public static explicit operator UserSurrogate(MemberShipUser other)
        {
            return  new UserSurrogate() { Name = other.Name };
        }
    }

    class User
    {
        public string Name { get; set; }
    }

    class MemberShipUser
    {
        public string Name { get; set; }   
    }


    }
}