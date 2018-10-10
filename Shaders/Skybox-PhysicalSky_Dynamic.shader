Shader "Custom/Skybox-PhysicalSky_Dynamic" 
{
	Properties 
	{
	}
	SubShader
	{
		Tags { "Queue"="Background" "RenderType"="Background" "PreviewType"="Skybox" }
		Cull Off ZWrite Off

		Pass
		{
			CGPROGRAM
			#pragma target 5.0
			#pragma vertex vert
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "Lighting.cginc"
			#include "PhysicalSkyCG.cginc"

			uniform float3 camera;
			uniform float3 sun_size;
			uniform float sky_exposure;
			uniform float sun_brightness;

			uniform samplerCUBE star_cubemap;
			uniform float star_brightness;
			uniform float4x4 star_rotation;

			#if defined(UNITY_COLORSPACE_GAMMA)
				#define LINEAR_TO_OUTPUT(color) sqrt(color)
			#else
				#define LINEAR_TO_OUTPUT(color) color
			#endif

			struct appdata_t
			{
				float4 vertex : POSITION;
				UNITY_VERTEX_INPUT_INSTANCE_ID
			};

			struct v2f
			{
				float4 pos : SV_POSITION;
				half3 view : TEXCOORD0;

				UNITY_VERTEX_OUTPUT_STEREO
			};

			v2f vert (appdata_t v)
			{
				v2f OUT;
				UNITY_SETUP_INSTANCE_ID(v);
				UNITY_INITIALIZE_VERTEX_OUTPUT_STEREO(OUT);
				OUT.pos = UnityObjectToClipPos(v.vertex);
				OUT.view = mul((float3x3)unity_ObjectToWorld, v.vertex.xyz);
				return OUT;
			}

			half4 frag (v2f IN) : SV_Target
			{
				float3 sun_direction = _WorldSpaceLightPos0.xyz;

				float3 view_ray = normalize(IN.view.xyz);
				float shadow_length = 0.0;

				float3 transmittance = 0;
				float3 luminance = GetSkyLuminance(camera, view_ray, shadow_length, sun_direction, transmittance);
				// HACK: No other real way of telling if we're being rendered as part of a reflection capture (in which we shouldn't be drawing the sun)
				bool reflection_capture = any(_LightColor0.xyz == half3(0, 0, 0));
				if (!reflection_capture && dot(sun_direction, view_ray) > sun_size.y)
				{
					luminance += transmittance * GetSolarLuminance();
				}

				luminance *= sun_brightness;

				luminance += texCUBE(star_cubemap, mul(view_ray, star_rotation)).rgb * star_brightness * transmittance;

				half3 result = LINEAR_TO_OUTPUT(luminance * sky_exposure);

				return half4(result, 1.0f);
			}
			ENDCG
		}

	}

	FallBack Off
}
