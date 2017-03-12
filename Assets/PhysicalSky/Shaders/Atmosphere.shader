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
				return float4(normalize(input.view.rgb) * 0.5 + 0.5, 0);
			}
			ENDCG
		}
		
	}
	Fallback Off
}
