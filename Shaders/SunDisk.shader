Shader "Hidden/Precompute"
{
	Properties {}
	SubShader
	{
		Cull Off ZWrite Off ZTest Always Fog{ Mode Off }

		// Pass 0, Compute Transmittance
		Pass
		{
			CGPROGRAM
			#pragma target 5.0
			#pragma vertex vert_img
			#pragma fragment frag
			
			#include "UnityCG.cginc"
			#include "PhysicalSkyCG.cginc"

			uniform float3 camera;
			uniform float3 sun_size;
			uniform float3 sun_direction;
			uniform float4x4 sun_rotation_matrix;

			fixed4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				// View rays, scaled to the size of the sun
				float3 untransformedviewray = normalize(float3(i.uv * 2 - 1, -(1 / tan(sun_size.x))));

				// Rotate the view rays to face towards the sun
				float3 view_ray = mul(sun_rotation_matrix, float4(untransformedviewray, 0)).xyz;
				float shadow_length = 0.0;

				float3 transmittance;
				GetSkyLuminance(camera, view_ray, shadow_length, sun_direction, transmittance);

				if (dot(sun_direction, view_ray) > sun_size.y)
				{
					return fixed4(transmittance * GetSolarRadiance(params), 1);
				}
				return 0;
			}
			ENDCG
		}
	}
	Fallback Off
}
