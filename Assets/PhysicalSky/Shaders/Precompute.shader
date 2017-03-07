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
			#include "PhysicalSkyCommon.cginc"
			#include "AtmosphereUniforms.cginc"

			fixed4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();
				float3 result = ComputeTransmittanceToTopAtmosphereBoundaryTexture(params, i.pos.xy);
				return fixed4(result.rgb, 1.0);
			}
			ENDCG
		}
		// Pass 1, Compute Direct Irradiance
		Pass
		{
			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"
			#include "AtmosphereUniforms.cginc"

			RWTexture2D<float4> delta_irradiance;

			uniform sampler2D transmittance_texture;

			float4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				float3 result = ComputeDirectIrradianceTexture(params, transmittance_texture, i.pos.xy);

				delta_irradiance[int2(i.pos.xy)] = float4(result, 0);
				return float4(0, 0, 0, 0); // Initializing irradianceLUT
			}
			ENDCG
		}
		// Pass 2, Compute Single Scattering
		Pass
		{
			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"
			#include "AtmosphereUniforms.cginc"

			RWTexture3D<float4> delta_rayleigh;
			RWTexture3D<float4> delta_mie;
			RWTexture3D<float4> scattering;

			uniform sampler2D transmittance_texture;

			fixed4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				int zz;
				for (zz = 0; zz < SCATTERING_TEXTURE_DEPTH; zz++)
				{
					float3 uvw = float3(i.pos.xy, zz + 0.5);
					float3 outRayleigh;
					float3 outMie;

					ComputeSingleScatteringTexture(params, transmittance_texture, uvw, outRayleigh, outMie);					

					int3 index = int3(i.pos.xy, zz);
					delta_rayleigh[index] = float4(outRayleigh, 1);
					delta_mie[index] = float4(outMie, 1);
					scattering[index] = float4(delta_rayleigh[index].rgb, delta_mie[index].r);
				}
				discard;
				return 0;
			}
			ENDCG
		}
		// Pass 3, Compute Scattering Density
		Pass
		{
			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"
			#include "AtmosphereUniforms.cginc"

			RWTexture3D<float4> scattering_density;

			uniform sampler2D transmittance_texture;
			uniform sampler3D single_rayleigh_scattering_texture;
			uniform sampler3D single_mie_scattering_texture;
			uniform sampler3D multiple_scattering_texture;
			uniform sampler2D irradiance_texture;
			uniform int scattering_order;

			fixed4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				int zz;
				for (zz = 0; zz < SCATTERING_TEXTURE_DEPTH; zz++)
				{
					float3 uvw = float3(i.pos.xy, zz + 0.5);

					float3 scatteringDensityResult = ComputeScatteringDensityTexture(
						params,
						transmittance_texture,
						single_rayleigh_scattering_texture,
						single_mie_scattering_texture,
						multiple_scattering_texture,
						irradiance_texture,
						uvw,
						scattering_order);

					// For debugging
					//scatteringDensityResult = float3(scattering_order == 2 ? 1 : 0, scattering_order == 3 ? 1 : 0, scattering_order == 4 ? 1 : 0);

					int3 index = int3(i.pos.xy, zz);
					scattering_density[index] = float4(scatteringDensityResult, 1);
				}
				discard;
				return 0;
			}
			ENDCG
		}
		// Pass 4, Compute Indirect Irradiance
		Pass
		{
			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"
			#include "AtmosphereUniforms.cginc"

			RWTexture3D<float4> volumeTarget;

			fixed4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				int zz;
				for (zz = 0; zz < SCATTERING_TEXTURE_DEPTH; zz++)
				{
					float3 pixel = float3(i.pos.xy, zz);

					float4 result = 0; // TODO

					volumeTarget[int3(i.pos.xy, zz)] = result;
				}
				discard;
				return 0;
			}
			ENDCG
		}
		// Pass 5, Compute Multiple Scattering
		Pass
		{
			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"
			#include "AtmosphereUniforms.cginc"

			RWTexture3D<float4> volumeTarget;

			fixed4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				int zz;
				for (zz = 0; zz < SCATTERING_TEXTURE_DEPTH; zz++)
				{
					float3 pixel = float3(i.pos.xy, zz);

					float4 result = 0; // TODO

					volumeTarget[int3(i.pos.xy, zz)] = result;
				}
				discard;
				return 0;
			}
			ENDCG
		}

	}
	Fallback Off
}
