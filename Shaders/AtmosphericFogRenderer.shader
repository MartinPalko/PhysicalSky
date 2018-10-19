Shader "Hidden/AtmosphericFogRenderer"
{
    HLSLINCLUDE

	#define HLSL_SHADER
    #include "PostProcessing/Shaders/StdLib.hlsl"
	#include "PhysicalSkyCG.cginc"

	#if defined(UNITY_COLORSPACE_GAMMA)
	#define LINEAR_TO_OUTPUT(color) sqrt(color)
	#else
	#define LINEAR_TO_OUTPUT(color) color
	#endif

	uniform float3 camera;
	uniform float3 sun_direction;
	uniform float sun_brightness;
	uniform float sky_exposure;
	uniform float4x4 _ClipToWorld;

	TEXTURE2D_SAMPLER2D(_CameraDepthTexture, sampler_CameraDepthTexture);
    TEXTURE2D_SAMPLER2D(_MainTex, sampler_MainTex);

	struct Attributes
	{
		float3 vertex : POSITION;
	};

	struct Varyings
	{
		float4 vertex : SV_POSITION;
		float2 texcoord : TEXCOORD0;
		float2 texcoordStereo : TEXCOORD1;
		float3 worldDirection : TEXCOORD2;
	};

	Varyings Vert(Attributes v)
	{
		Varyings o;
		o.vertex = float4(v.vertex.xy, 0.0, 1.0);
		o.texcoord = TransformTriangleVertexToUV(v.vertex.xy);

#if UNITY_UV_STARTS_AT_TOP
		o.texcoord = o.texcoord * float2(1.0, -1.0) + float2(0.0, 1.0);
#endif
		o.texcoordStereo = TransformStereoScreenSpaceTex(o.texcoord, 1.0);
		o.worldDirection = mul(_ClipToWorld, o.vertex) - _WorldSpaceCameraPos;
		return o;
	}

    float4 Frag(Varyings i) : SV_Target
    {
        float4 sceneColor = SAMPLE_TEXTURE2D(_MainTex, sampler_MainTex, i.texcoord);

		float rawDepth = SAMPLE_DEPTH_TEXTURE(_CameraDepthTexture, sampler_CameraDepthTexture, i.texcoord);
		float depth = LinearEyeDepth(rawDepth);
		float3 wDir = i.worldDirection * depth;

		// Don't double-scatter on skybox
		[branch]
#if defined(UNITY_REVERSED_Z)
		if(rawDepth < 0.0000001f)
#else
		if(rawDepth > 0.9999999f)
#endif
#ifndef USE_FINAL_BLEND
			return sceneColor;
#else
			clip(-1);
#endif

		const float3 shadow_length = 0;

		float3 transmittance;
		float3 in_scatter = GetSkyLuminanceToPoint(camera, normalize(i.worldDirection), length(wDir / 1000), shadow_length, sun_direction, transmittance);

		in_scatter *= sun_brightness;

		in_scatter = LINEAR_TO_OUTPUT(in_scatter * sky_exposure);

		return float4(sceneColor.rgb * transmittance + in_scatter, sceneColor.a);
    }

    ENDHLSL

    SubShader
    {
        Cull Off ZWrite Off ZTest Always

        Pass
        {
            HLSLPROGRAM

                #pragma vertex Vert
                #pragma fragment Frag

            ENDHLSL
        }
    }

}