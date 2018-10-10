

using UnityEditor;
using UnityEditor.AnimatedValues;
using UnityEngine;

[CustomEditor(typeof(PhysicalSky.AtmosphereModel))]
public class AtmosphereModelEditor : Editor
{
    // Animated values
    protected AnimBool m_IsShowingAdvancedRegion;

    // Content
    private GUIContent m_AdvancedModeLabel;

    // Properties
    protected SerializedProperty m_Script;
    protected SerializedProperty m_scatteringOrders;
    protected SerializedProperty m_luminance;
    protected SerializedProperty m_useOzone;
    protected SerializedProperty m_useConstantSolarSpectrum;
    protected SerializedProperty m_constantSolarIrradiance;
    protected SerializedProperty m_sunAngularRadius;
    protected SerializedProperty m_planetaryRadius;
    protected SerializedProperty m_atmosphereThickness;
    protected SerializedProperty m_rayleigh;
    protected SerializedProperty m_rayleighScaleHeight;
    protected SerializedProperty m_mieScaleHeight;
    protected SerializedProperty m_mieAngstromAlpha;
    protected SerializedProperty m_mieAngstromBeta;
    protected SerializedProperty m_mieSingleScatteringAlbedo;
    protected SerializedProperty m_miePhaseFunctionG;
    protected SerializedProperty m_groundAlbedo;
    protected SerializedProperty m_maxSunZenithAngle;
    protected SerializedProperty m_PrecomputeShader;

    private bool isAdvancedModeExpanded
    {
        get { return m_Script.isExpanded; }
        set { m_Script.isExpanded = value; }
    }

    /// <summary>
    /// This function is called when the object is loaded.
    /// We use it to init all our properties
    /// </summary>
    protected virtual void OnEnable()
    {
        m_Script = serializedObject.FindProperty("m_Script");
        m_scatteringOrders = serializedObject.FindProperty("m_scatteringOrders");
        m_luminance = serializedObject.FindProperty("m_luminance");
        m_useOzone = serializedObject.FindProperty("m_useOzone");
        m_useConstantSolarSpectrum = serializedObject.FindProperty("m_useConstantSolarSpectrum");
        m_constantSolarIrradiance = serializedObject.FindProperty("m_constantSolarIrradiance");
        m_sunAngularRadius = serializedObject.FindProperty("m_sunAngularRadius");
        m_planetaryRadius = serializedObject.FindProperty("m_planetaryRadius");
        m_atmosphereThickness = serializedObject.FindProperty("m_atmosphereThickness");
        m_rayleigh = serializedObject.FindProperty("m_rayleigh");
        m_rayleighScaleHeight = serializedObject.FindProperty("m_rayleighScaleHeight");
        m_mieScaleHeight = serializedObject.FindProperty("m_mieScaleHeight");
        m_mieAngstromAlpha = serializedObject.FindProperty("m_mieAngstromAlpha");
        m_mieAngstromBeta = serializedObject.FindProperty("m_mieAngstromBeta");
        m_mieSingleScatteringAlbedo = serializedObject.FindProperty("m_mieSingleScatteringAlbedo");
        m_miePhaseFunctionG = serializedObject.FindProperty("m_miePhaseFunctionG");
        m_groundAlbedo = serializedObject.FindProperty("m_groundAlbedo");
        m_maxSunZenithAngle = serializedObject.FindProperty("m_maxSunZenithAngle");
        m_PrecomputeShader = serializedObject.FindProperty("m_PrecomputeShader");

        m_IsShowingAdvancedRegion = new AnimBool(isAdvancedModeExpanded, Repaint);

        m_AdvancedModeLabel = new GUIContent("Advanced Settings", "Exposes all the settings that are used internally.");
    }

    protected virtual void OnDisable()
    {
        m_IsShowingAdvancedRegion.valueChanged.RemoveListener(Repaint);
        m_IsShowingAdvancedRegion = null;
    }

    /// <summary>
    /// Inside this function you can add your own custom GUI
    /// for the inspector of a specific object class.
    /// </summary>
    public override void OnInspectorGUI()
    {
        EditorGUI.BeginChangeCheck();
        {
            EditorGUILayout.PropertyField(m_Script);
            isAdvancedModeExpanded = EditorGUILayout.Foldout(isAdvancedModeExpanded, m_AdvancedModeLabel);
            if (EditorGUILayout.BeginFadeGroup(m_IsShowingAdvancedRegion.faded))
            {
                EditorGUI.indentLevel++;
                EditorGUILayout.PropertyField(property: m_scatteringOrders);
                EditorGUILayout.PropertyField(property: m_luminance);
                EditorGUILayout.PropertyField(property: m_useOzone);
                EditorGUILayout.PropertyField(property: m_useConstantSolarSpectrum);
                EditorGUILayout.PropertyField(property: m_constantSolarIrradiance);
                EditorGUILayout.PropertyField(property: m_sunAngularRadius);
                EditorGUILayout.PropertyField(property: m_planetaryRadius);
                EditorGUILayout.PropertyField(property: m_atmosphereThickness);
                EditorGUILayout.PropertyField(property: m_rayleigh);
                EditorGUILayout.PropertyField(property: m_rayleighScaleHeight);
                EditorGUILayout.PropertyField(property: m_mieScaleHeight);
                EditorGUILayout.PropertyField(property: m_mieAngstromAlpha);
                EditorGUILayout.PropertyField(property: m_mieAngstromBeta);
                EditorGUILayout.PropertyField(property: m_mieSingleScatteringAlbedo);
                EditorGUILayout.PropertyField(property: m_miePhaseFunctionG);
                EditorGUILayout.PropertyField(property: m_groundAlbedo);
                EditorGUILayout.PropertyField(property: m_maxSunZenithAngle);
                EditorGUILayout.PropertyField(property: m_PrecomputeShader);
                EditorGUI.indentLevel--;
            }
            EditorGUILayout.EndFadeGroup();
        }
        if (EditorGUI.EndChangeCheck())
        {
            serializedObject.ApplyModifiedProperties();

            PhysicalSky.AtmosphereModel model = serializedObject.targetObject as PhysicalSky.AtmosphereModel;
            if (model != null)
                model.Compute();
        }

        m_IsShowingAdvancedRegion.target = isAdvancedModeExpanded;
    }
}
