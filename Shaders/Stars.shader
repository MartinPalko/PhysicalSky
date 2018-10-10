Shader "Stars" 
{
	Properties 
	{
	}
	SubShader
	{
		Tags { "Queue"="Transparent" "RenderType"="Background" }
		Cull Off ZWrite Off
		Blend One One

		Pass
		{
			CGPROGRAM
			#pragma target 5.0
			#pragma vertex vert
			#pragma geometry geom
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "Lighting.cginc"
			#include "PhysicalSkyCommon.cginc"

			uniform sampler2D transmittance_texture;

			uniform float3 camera;

			uniform float star_intensity_multiplier;
			uniform float star_intensity_power;
			uniform float planet_size;

			#if defined(UNITY_COLORSPACE_GAMMA)
			#define LINEAR_TO_OUTPUT(color) sqrt(color)
			#else
			#define LINEAR_TO_OUTPUT(color) color
			#endif

			struct appdata_t
			{
				float4 vertex : POSITION;
				float4 color : COLOR;
				UNITY_VERTEX_INPUT_INSTANCE_ID
			};

			struct v2f
			{
				float4 pos : SV_POSITION;
				half2 uv : TEXCOORD0;
				float4 color : COLOR;
				UNITY_VERTEX_OUTPUT_STEREO
			};

			v2f vert (appdata_t v)
			{
				v2f OUT;
				UNITY_SETUP_INSTANCE_ID(v);
				UNITY_INITIALIZE_VERTEX_OUTPUT_STEREO(OUT);
				OUT.uv = 0;

				OUT.pos = v.vertex;
				OUT.pos.w = 0; // No translation for model and view transformations; always be centered around 0
				OUT.pos = mul(UNITY_MATRIX_M, OUT.pos); // To world space

				const float3 worldDirection = normalize(OUT.pos).xyz; // Save for later.

				OUT.pos = mul(UNITY_MATRIX_V, OUT.pos); // To view space
				OUT.pos.w = 1; // Translation is needed again for projection
				OUT.pos = mul(UNITY_MATRIX_P, OUT.pos); // Project

				AtmosphereParameters params = GetAtmosphereParameters();

				// Set color to 0 if we're below the horizon.
				float angle = acos(dot(float3(0,-1,0), worldDirection)) * 2;
				OUT.color = angle > planet_size ? v.color : 0;

				if (OUT.color.a > 0)
				{
					float r = length(camera);
					float rmu = dot(camera, worldDirection);
					float mu = rmu / r;
					float3 atmosphereTransmittance = GetTransmittanceToTopAtmosphereBoundary(params, transmittance_texture, r, mu);
					OUT.color.rgb *= atmosphereTransmittance;
				}

				return OUT;
			}

			[maxvertexcount(4 * 3)] // 1 quad per vertex
			void geom(triangle v2f input[3], inout TriangleStream<v2f> OutputStream)
			{
				v2f newVert = (v2f)0;

				float2 quadSize = 2.8f / _ScreenParams.xy;

				// Far Z value, to make sure starts are always at the end of the frustum
				float farZ = 0;
				// Put z value just on the far plane
				if (UNITY_NEAR_CLIP_VALUE == -1)
					farZ = -1e-7; // Just below 0
				else if (UNITY_NEAR_CLIP_VALUE == 1)
					farZ = 1e-7; // Just above 0
				else
					farZ = 1 - 1e-7; // Just below 1

				[unroll]
				for(int i = 0; i < 3; i++)
				{
					if (input[i].color.a > 0)
					{
						newVert.color = input[i].color;

						input[i].pos.xyz /= input[i].pos.w;

						// Top Left
						newVert.pos = input[i].pos + float4(-quadSize.x, quadSize.y, 0, 0);
						newVert.pos.xy = newVert.pos * newVert.pos.w;
						newVert.pos.z = farZ;
						newVert.uv = float2(-1, 1);
						OutputStream.Append(newVert);

						// Bottom Left
						newVert.pos = input[i].pos + float4(-quadSize.x, -quadSize.y, 0, 0);
						newVert.pos.xy = newVert.pos * newVert.pos.w;
						newVert.pos.z = farZ;
						newVert.uv = float2(-1, -1);
						OutputStream.Append(newVert);

						// Top Right
						newVert.pos = input[i].pos + float4(quadSize.x, quadSize.y, 0, 0);
						newVert.pos.xy = newVert.pos * newVert.pos.w;
						newVert.pos.z = farZ;
						newVert.uv = float2(1, 1);
						OutputStream.Append(newVert);

						// Bottom Right
						newVert.pos = input[i].pos + float4(quadSize.x, -quadSize.y, 0, 0);
						newVert.pos.xy = newVert.pos * newVert.pos.w;
						newVert.pos.z = farZ;
						newVert.uv = float2(1, -1);
						OutputStream.Append(newVert);

						OutputStream.RestartStrip();
					}
				}
			}

			half4 frag (v2f IN) : SV_Target
			{
				half falloff = saturate(1 - length(IN.uv));
				half3 color = IN.color.rgb;
				half3 intensity = pow(IN.color.a, star_intensity_power) * star_intensity_multiplier;

				half3 result = color * intensity * falloff;
				return half4(LINEAR_TO_OUTPUT(result), 1);
			}
			ENDCG	
		}

	}

	FallBack Off
}