Shader "Custom/Atmosphere" 
{
	Properties 
	{
	}
	SubShader 
	{
		Cull Off ZWrite Off ZTest LEqual Fog{ Mode Off }

		Pass
		{
			CGPROGRAM
			#pragma target 3.0
			#pragma vertex vert
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"
			#include "AtmosphereUniforms.cginc"

			sampler2D _MainTex;
			uniform float4 _MainTex_TexelSize;
			uniform sampler2D transmittance_texture;
			uniform sampler3D scattering_texture;
			uniform sampler2D irradiance_texture;

			uniform float3 camera;
			uniform float3 sun_direction;
			uniform float3 sun_radiance;
			uniform float3 sun_size;

			struct appdata_t
			{
				float4 pos : POSITION;
			};

			struct v2f
			{
				float4 pos : SV_POSITION;
				float4 view : TEXCOORD0;
			};

			float4 real_cam_pos;

			v2f vert(appdata_t input)
			{
				v2f output;
				output.view = mul(unity_ObjectToWorld, input.pos);
				output.pos = input.pos;

				// Put z value just on the far plane
				if (UNITY_NEAR_CLIP_VALUE == -1)
					output.pos.z = -1e-7; // Just below 0
				else if (UNITY_NEAR_CLIP_VALUE == 1)
					output.pos.z = 1e-7; // Just above 0
				else
					output.pos.z = 1 - 1e-7; // Just below 1

				// Compensate for flipped projection matrix.
				if (_ProjectionParams.x < 0)
					output.pos.y *= -1;

				return output;
			}

			half4 frag(v2f input) : COLOR
			{
				// For Debugging
				//return float4(normalize(input.view.rgb) * 0.5 + 0.5, 0);

				AtmosphereParameters params = GetAtmosphereParameters();

				float3 view_ray = normalize(input.view.xyz);
				float lightshaft_fadein_hack = smoothstep(0.02, 0.04, dot(camera, sun_direction));
				float shadow_length = 0.0;


				float3 transmittance;
				float3 radiance = GetSkyRadiance(params, transmittance_texture, scattering_texture, scattering_texture, camera, view_ray, shadow_length, sun_direction, transmittance);

				if (dot(sun_direction, view_ray) > sun_size.y)
				{
					radiance += transmittance * sun_radiance;
				}

				float exposure = 10;
				float white_point = 0.4;

				return float4(pow(1 - exp(-radiance / white_point * exposure), (1.0 / 2.2)), 1);
			}
			ENDCG
		}
		
	}
	Fallback Off
}
