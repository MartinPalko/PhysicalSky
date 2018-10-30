using System.Linq;
using UnityEngine;
using UnityEditor;
using System.IO;

namespace PhysicalSky
{
    static class ShaderIncludes
    {
#if UNITY_2018_1_OR_NEWER
        [ShaderIncludePath]
#endif
        public static string[] GetPaths()
        {
            return new string[] {
                Path.Combine(Application.dataPath, "PhysicalSky/Shaders/"),
                };
        }
    }
}
