using UnityEngine;
using UnityEditor;
using System;

namespace PhysicalSky
{
    [CustomEditor(typeof(AtmosphereModel))]
    public class AtmosphereModelEditor : Editor
    {
        // Toggles are declared static, so they don't keep collapsing between entering/exiting the editor.
        static bool advancedUnfolded = false;

        public override void OnInspectorGUI()
        {
            AtmosphereModel atmosphereModel = serializedObject.targetObject as AtmosphereModel;

            if (!atmosphereModel)
                return;

            advancedUnfolded = EditorGUILayout.Foldout(advancedUnfolded, "Advanced", true);
            if (advancedUnfolded)
            {
                atmosphereModel.ConstantSolarIrradiance = Mathf.Max(0.001f, EditorGUILayout.FloatField("Solar Irradiance", atmosphereModel.ConstantSolarIrradiance));
                atmosphereModel.SunAngularRadius = Mathf.Max(0.001f, EditorGUILayout.FloatField("Sun Radius (deg)", atmosphereModel.SunAngularRadius * Mathf.PI * 2.0f)) / (Mathf.PI * 2.0f);
                atmosphereModel.PlanetaryRadius = Mathf.Max(1.0f, EditorGUILayout.FloatField("Planetary Radius (km)", atmosphereModel.PlanetaryRadius / 1000) * 1000);
                atmosphereModel.AtmosphereThickness = Mathf.Max(1.0f, EditorGUILayout.FloatField("Atmosphere Thickness (km)", atmosphereModel.AtmosphereThickness / 1000) * 1000);
                atmosphereModel.Rayleigh = Mathf.Max(1.0e-22f, EditorGUILayout.FloatField("Rayleigh", atmosphereModel.Rayleigh));
                atmosphereModel.RayleighScaleHeight = Mathf.Max(0.0f, EditorGUILayout.FloatField("RayleighScaleHeight", atmosphereModel.RayleighScaleHeight));
                atmosphereModel.MieScaleHeight = Mathf.Max(0.0f, EditorGUILayout.FloatField("MieScaleHeight", atmosphereModel.MieScaleHeight));
                atmosphereModel.MieAngstromAlpha = Mathf.Max(0.0f, EditorGUILayout.FloatField("MieAngstromAlpha", atmosphereModel.MieAngstromAlpha));
                atmosphereModel.MieAngstromBeta = Mathf.Max(0.0f, EditorGUILayout.FloatField("MieAngstromBeta", atmosphereModel.MieAngstromBeta));
                atmosphereModel.MieSingleScatteringAlbedo = Mathf.Max(0.0f, EditorGUILayout.FloatField("MieSingleScatteringAlbedo", atmosphereModel.MieSingleScatteringAlbedo));
                atmosphereModel.MiePhaseFunctionG = Mathf.Max(0.0f, EditorGUILayout.FloatField("MiePhaseFunctionG", atmosphereModel.MiePhaseFunctionG));
                atmosphereModel.GroundAlbedo = Mathf.Max(0.0f, EditorGUILayout.FloatField("GroundAlbedo", atmosphereModel.GroundAlbedo));
                atmosphereModel.MaxSunZenithAngle = Mathf.Max(0.0f, EditorGUILayout.FloatField("MaxSunZenithAngle", atmosphereModel.MaxSunZenithAngle));
                atmosphereModel.PrecomputeShader = EditorGUILayout.ObjectField("Precompute Shader", atmosphereModel.PrecomputeShader as Shader, typeof(Shader), false) as Shader;
            }

            EditorUtility.SetDirty(atmosphereModel);

        }
    }
}