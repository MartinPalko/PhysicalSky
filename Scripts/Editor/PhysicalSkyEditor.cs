

using UnityEditor;
using UnityEditor.AnimatedValues;
using UnityEngine;

[CanEditMultipleObjects]
[CustomEditor(typeof(PhysicalSky.PhysicalSky))]
public class PhysicalSkyEditor : Editor
{
    public class Styles
    {
        public GUIStyle regionHeader;

        public Styles()
        {
            regionHeader = new GUIStyle(EditorStyles.foldout);
            regionHeader.fixedHeight = EditorGUIUtility.singleLineHeight;
            regionHeader.stretchWidth = true;
            regionHeader.fontStyle = FontStyle.Bold;
        }
    }


    // Animations
    private AnimBool m_AnimAdvancedRegion;
    private AnimBool m_AnimSunRegion;
    private AnimBool m_AnimStarRegion;

    
    // Styles
    private Styles m_Style;

    // Content
    private GUIContent m_SunRegionLabel;
    private GUIContent m_StarRegionLabel;
    private GUIContent m_AdvancedRegionLabel;

    // Properties
    protected SerializedProperty m_Script;
    protected SerializedProperty m_SunLightBrightnessMultiplier;
    protected SerializedProperty m_SunBrightnessMultiplier;
    protected SerializedProperty m_SunSaturationMultiplier;
    protected SerializedProperty m_SkyExposure;
    protected SerializedProperty m_Altitude;
    protected SerializedProperty m_Atmosphere;
    protected SerializedProperty m_StarBrightnessMultiplier;
    protected SerializedProperty m_StarBrightnessPower;
    protected SerializedProperty m_StarMap;
    protected SerializedProperty m_SkyShader;
    protected SerializedProperty m_SunRadianceShader;
    protected SerializedProperty m_StarMeshShader;

    /// <summary>
    /// This function is called when the object is loaded.
    /// We use it to init all our properties
    /// </summary>
    protected virtual void OnEnable()
    {
        m_SunLightBrightnessMultiplier = serializedObject.FindProperty("m_SunLightBrightnessMultiplier");
        m_SunBrightnessMultiplier = serializedObject.FindProperty("m_SunBrightnessMultiplier");
        m_SunSaturationMultiplier = serializedObject.FindProperty("m_SunSaturationMultiplier");
        m_SkyExposure = serializedObject.FindProperty("m_SkyExposure");
        m_Altitude = serializedObject.FindProperty("m_Altitude");
        m_Atmosphere = serializedObject.FindProperty("m_Atmosphere");
        m_StarBrightnessMultiplier = serializedObject.FindProperty("m_StarBrightnessMultiplier");
        m_StarBrightnessPower = serializedObject.FindProperty("m_StarBrightnessPower");
        m_StarMap = serializedObject.FindProperty("m_StarMap");
        m_SkyShader = serializedObject.FindProperty("m_SkyShader");
        m_SunRadianceShader = serializedObject.FindProperty("m_SunRadianceShader");
        m_StarMeshShader = serializedObject.FindProperty("m_StarMeshShader");

        // Labels
        m_AdvancedRegionLabel = new GUIContent("Advanced Settings", "This contains a few more advanced settings not required by most users");
        m_SunRegionLabel = new GUIContent("Sun Settings", "These are all the settings for the sun.");
        m_StarRegionLabel = new GUIContent("Star Settings", "These are all the settings for the stars.");

        // We use the isExpanded from the m_StarMeshShader but we could really use any.
        // This is saved by Unity so it survives compiles.
        m_AnimAdvancedRegion = new AnimBool(isAdvancedRegionExpanded, Repaint);
        m_AnimStarRegion = new AnimBool(isStarRegionExpanded, Repaint);
        m_AnimSunRegion = new AnimBool(IsSunRegionExpanded, Repaint);
    }

    protected void OnDisable()
    {
        m_AnimAdvancedRegion.valueChanged.RemoveListener(Repaint);
        m_AnimStarRegion.valueChanged.RemoveListener(Repaint);
        m_AnimSunRegion.valueChanged.RemoveListener(Repaint);
    }

    /// <summary>
    /// We use isExpaned from any property in this class to see
    /// if our foldout is open. This is serialized by Unity
    /// behind the scenes so it survives a compile.
    /// </summary>
    private bool isAdvancedRegionExpanded
    {
        get { return m_StarMeshShader.isExpanded; }
        set { m_StarMeshShader.isExpanded = value; }
    }

    private bool IsSunRegionExpanded
    {
        get { return m_SunLightBrightnessMultiplier.isExpanded; }
        set { m_SunLightBrightnessMultiplier.isExpanded = value; }
    }

    private bool isStarRegionExpanded
    {
        get { return m_StarMap.isExpanded; }
        set { m_StarMap.isExpanded = value; }
    }

    /// <summary>
    /// Inside this function you can add your own custom GUI
    /// for the inspector of a specific object class.
    /// </summary>
    public override void OnInspectorGUI()
    {
        if(m_Style == null)
        {
            // Styles can only be created during an OnGUI callback.
            m_Style = new Styles();
        }

        EditorGUI.BeginChangeCheck();
        {
            EditorGUILayout.PropertyField(property: m_Altitude);
            EditorGUILayout.PropertyField(property: m_Atmosphere);
            EditorGUILayout.PropertyField(property: m_SkyExposure);

            IsSunRegionExpanded = EditorGUILayout.Foldout(IsSunRegionExpanded, m_SunRegionLabel);
            if (EditorGUILayout.BeginFadeGroup(m_AnimSunRegion.faded))
            {
                EditorGUI.indentLevel++;
                EditorGUILayout.PropertyField(property: m_SunLightBrightnessMultiplier);
                EditorGUILayout.PropertyField(property: m_SunBrightnessMultiplier);
                EditorGUILayout.PropertyField(property: m_SunSaturationMultiplier);
                EditorGUI.indentLevel--;
            }
            EditorGUILayout.EndFadeGroup();

            isStarRegionExpanded = EditorGUILayout.Foldout(isStarRegionExpanded, m_StarRegionLabel);
            if (EditorGUILayout.BeginFadeGroup(m_AnimStarRegion.faded))
            {
                EditorGUI.indentLevel++;
                EditorGUILayout.PropertyField(property: m_StarBrightnessMultiplier);
                EditorGUILayout.PropertyField(property: m_StarBrightnessPower);
                EditorGUILayout.PropertyField(property: m_StarMap);
                EditorGUI.indentLevel--;
            }
            EditorGUILayout.EndFadeGroup();

            isAdvancedRegionExpanded = EditorGUILayout.Foldout(isAdvancedRegionExpanded, m_AdvancedRegionLabel);
            if (EditorGUILayout.BeginFadeGroup(m_AnimAdvancedRegion.faded))
            {
                EditorGUI.indentLevel++;
                EditorGUILayout.PropertyField(property: m_SkyShader);
                EditorGUILayout.PropertyField(property: m_SunRadianceShader);
                EditorGUILayout.PropertyField(property: m_StarMeshShader);
                EditorGUI.indentLevel--;
            }
            EditorGUILayout.EndFadeGroup();
        }
        if (EditorGUI.EndChangeCheck())
        {
            serializedObject.ApplyModifiedProperties();
        }

        // Update our animation target
        m_AnimAdvancedRegion.target = isAdvancedRegionExpanded;
        m_AnimSunRegion.target = IsSunRegionExpanded;
        m_AnimStarRegion.target = isStarRegionExpanded;
    }
}
