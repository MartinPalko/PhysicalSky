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
			#pragma vertex vert
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "Lighting.cginc"
			#include "PhysicalSkyCommon.cginc"
			#include "AtmosphereUniforms.cginc"

			sampler2D _MainTex;
			uniform float4 _MainTex_TexelSize;
			uniform sampler2D transmittance_texture;
			uniform sampler3D scattering_texture;
			uniform sampler2D irradiance_texture;

			uniform float3 camera;
			uniform float3 sun_radiance;
			uniform float3 sun_size;

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

				AtmosphereParameters params = GetAtmosphereParameters();

				float3 view_ray = normalize(IN.view.xyz);
				float lightshaft_fadein_hack = smoothstep(0.02, 0.04, dot(camera, sun_direction));
				float shadow_length = 0.0;

				float3 transmittance;
				float3 radiance = GetSkyRadiance(params, transmittance_texture, scattering_texture, scattering_texture, camera, view_ray, shadow_length, sun_direction, transmittance);

				// HACK: No other real way of telling if we're being rendered as part of a reflection capture (in which we shouldn't be drawing the sun)
				bool reflection_capture = any(_LightColor0.xyz > half3(0, 0, 0));
				if (!reflection_capture && dot(sun_direction, view_ray) > sun_size.y)
				{
					radiance += transmittance * sun_radiance;
				}

				// TODO: Make into a uniform and set from script.
				float exposure = 20;
				float white_point = 0.4;

				return half4(radiance * exposure, 1.0f);
				//return half4(pow(1 - exp(-radiance / white_point * exposure), (1.0 / 2.2)), 1);

			}
			ENDCG
		}

	}

	FallBack Off
}
