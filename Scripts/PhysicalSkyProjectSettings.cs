using UnityEngine;

namespace PhysicalSky
{
    public sealed class PhysicalSkyProjectSettings : ScriptableObject
    {
        private const string ASSET_NAME = "PhysicalSkyProjectSettings";
        private const string SHADERS_FOLDER = "Assets/PhysicalSky/Shaders";
        private static PhysicalSkyProjectSettings instance = null;

        public static PhysicalSkyProjectSettings GetInstance()
        {
            if (instance)
                return instance;

            instance = Resources.Load<PhysicalSkyProjectSettings>(ASSET_NAME);

            if (instance)
                return instance;

#if UNITY_EDITOR
            const string savePath = "Assets/Resources/" + ASSET_NAME + ".asset";
            Debug.Log("Physical Sky Project Settings not found, creating at \"" + savePath + "\"");

            instance = CreateInstance<PhysicalSkyProjectSettings>();
            if (!UnityEditor.AssetDatabase.IsValidFolder("Assets/Resources/"))
                UnityEditor.AssetDatabase.CreateFolder("Assets", "Resources");
            UnityEditor.AssetDatabase.CreateAsset(instance, savePath);
            UnityEditor.AssetDatabase.SaveAssets();
            UnityEditor.AssetDatabase.Refresh();
            
            return instance;
#else
            throw new System.Exception("No PhysicalSkyProjectSettings have been included in the build's resources!");
#endif
        }

#if UNITY_EDITOR
        // Open Project Text Settings
        [UnityEditor.MenuItem("Edit/Project Settings/Physical Sky")]
        public static void SelectProjectSettings()
        {
            UnityEditor.Selection.activeObject = GetInstance();
        }
#endif

        // Make constructor private
        private PhysicalSkyProjectSettings() { }

        [SerializeField]
        private Shader m_prerenderBlitShader = null;
        public static Shader PrerenderBlitShader
        {
            get { return GetInstance().m_prerenderBlitShader; }
            set { GetInstance().m_prerenderBlitShader = value; }
        }

        [SerializeField]
        private ComputeShader m_prerenderComputeShader = null;
        public static ComputeShader PrerenderComputeShader
        {
            get { return GetInstance().m_prerenderComputeShader; }
            set { GetInstance().m_prerenderComputeShader = value; }
        }

    }
}
