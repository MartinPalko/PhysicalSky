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

			float4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();
				float3 result = ComputeTransmittanceToTopAtmosphereBoundaryTexture(params, i.pos.xy);
				return float4(result.rgb, 1.0);
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

			uniform sampler2D transmittance_texture;

			struct f2a
			{
				float4 delta_irradiance : COLOR0;
				float4 irradiance_lut : COLOR1;
			};

			f2a frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				float3 result = ComputeDirectIrradianceTexture(params, transmittance_texture, i.pos.xy);

				f2a OUT;
				OUT.delta_irradiance = float4(result, 0);
				OUT.irradiance_lut =  float4(0, 0, 0, 0); // Initializing irradianceLUT
				return OUT;
			}
			ENDCG
		}
		// Pass 2, Compute Single Scattering
		Pass
		{
			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img3d
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "Blit3d.cginc"
			#include "PhysicalSkyCommon.cginc"

			uniform sampler2D transmittance_texture;

			struct f2a
			{
				float4 delta_rayleigh : COLOR0;
				float4 delta_mie : COLOR1;
				float4 scattering : COLOR2;
			};

			f2a frag (v2f_img3d i)
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				float3 uvw = float3(i.pos.xy, i.uvw.z * SCATTERING_TEXTURE_DEPTH);
				float3 outRayleigh;
				float3 outMie;

				ComputeSingleScatteringTexture(params, transmittance_texture, uvw, outRayleigh, outMie);

				f2a OUT;
				OUT.delta_rayleigh = float4(outRayleigh, 1);
				OUT.delta_mie = float4(outMie, 1);
				OUT.scattering = float4(0, 0, 0, 0); // Initializing scatteringLUT
				return OUT;
			}
			ENDCG
		}
		// Pass 3, Combine Single Scattering
		Pass
		{
			Tags {"RenderType"="Transparent"}
			Blend One One // Additive Blend

			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img3d
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "Blit3d.cginc"
			#include "PhysicalSkyCommon.cginc"

			uniform sampler3D single_rayleigh_scattering_texture;
			uniform sampler3D single_mie_scattering_texture;

			uniform float _luminance_from_radiance[9];

			float4 frag (v2f_img3d i) : COLOR
			{
				float3x3 luminance_from_radiance = {_luminance_from_radiance};

				float3 rayleigh = tex3Dlod(single_rayleigh_scattering_texture, float4(i.uvw, 0)).rgb;
				float3 mie = tex3Dlod(single_mie_scattering_texture, float4(i.uvw, 0)).rgb;

				return float4(mul(luminance_from_radiance, rayleigh), mul(luminance_from_radiance, mie).r);
			}
			ENDCG
		}
		// Pass 4, Compute Scattering Density
		Pass
		{
			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img3d
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "Blit3d.cginc"
			#include "PhysicalSkyCommon.cginc"

			uniform sampler2D transmittance_texture;
			uniform sampler3D single_rayleigh_scattering_texture;
			uniform sampler3D single_mie_scattering_texture;
			uniform sampler3D multiple_scattering_texture;
			uniform sampler2D irradiance_texture;

			uniform int scattering_order;

			float4 frag (v2f_img3d i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				float3 uvw = float3(i.pos.xy, i.uvw.z * SCATTERING_TEXTURE_DEPTH);

				float3 scatteringDensityResult = ComputeScatteringDensityTexture(
					params,
					transmittance_texture,
					single_rayleigh_scattering_texture,
					single_mie_scattering_texture,
					multiple_scattering_texture,
					irradiance_texture,
					uvw,
					scattering_order);
					
				return float4(scatteringDensityResult, 1);
			}
			ENDCG
		}
		// Pass 5, Compute Indirect Irradiance
		Pass
		{
			CGPROGRAM
			
			#pragma target 5.0

			#pragma vertex vert_img
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"

			uniform sampler3D single_rayleigh_scattering_texture;
			uniform sampler3D single_mie_scattering_texture;
			uniform sampler3D multiple_scattering_texture;
			uniform sampler2D irradiance_texture;
			
			uniform int scattering_order;

			float4 frag (v2f_img i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				int2 index = i.pos.xy;
				float3 irradianceResult = ComputeIndirectIrradianceTexture(
					params,
					single_rayleigh_scattering_texture,
					single_mie_scattering_texture,
					multiple_scattering_texture,
					i.pos.xy,
					scattering_order);

				return float4(irradianceResult, 0);
			}
			ENDCG
		}
		// Pass 6, Accumulate Indirect Irradiance
		Pass
		{
			Tags {"RenderType"="Transparent"}
			Blend One One // Additive Blend

			CGPROGRAM
			
			#pragma target 5.0

			#pragma vertex vert_img
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"

			uniform sampler2D irradiance_texture;

			uniform float _luminance_from_radiance[9];

			float4 frag (v2f_img i) : COLOR
			{
				float3x3 luminance_from_radiance = {_luminance_from_radiance};

				return float4(mul(luminance_from_radiance, tex2Dlod(irradiance_texture, float4(i.uv, 0, 0)).rgb), 0);
			}
			ENDCG
		}
		// Pass 7, Compute Multiple Scattering
		Pass
		{
			CGPROGRAM
			
			#pragma target 5.0

			#pragma vertex vert_img3d
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "Blit3d.cginc"
			#include "PhysicalSkyCommon.cginc"

			uniform sampler2D transmittance_texture;
			uniform sampler3D scattering_density_texture;

			float4 frag (v2f_img3d i) : COLOR
			{
				AtmosphereParameters params = GetAtmosphereParameters();

				float3 uvw = float3(i.pos.xy, i.uvw.z * SCATTERING_TEXTURE_DEPTH);

				float nu;
				float3 deltaMultipleScatteringResult = ComputeMultipleScatteringTexture(
					params,
					transmittance_texture,
					scattering_density_texture,
					uvw,
					nu);

				return float4(deltaMultipleScatteringResult, nu);
			}
			ENDCG
		}
		// Pass 8, Accumulate Multiple Scattering
		Pass
		{
			Tags {"RenderType"="Transparent"}
			Blend One One // Additive Blend
			CGPROGRAM
			
			#pragma target 5.0

			#pragma vertex vert_img3d
			#pragma fragment frag

			#include "UnityCG.cginc"
			#include "Blit3d.cginc"
			#include "PhysicalSkyCommon.cginc"

			uniform sampler3D multiple_scattering_texture;

			uniform float _luminance_from_radiance[9];

			float4 frag (v2f_img3d i) : COLOR
			{
				float3x3 luminance_from_radiance = {_luminance_from_radiance};

				float4 texSample = tex3Dlod(multiple_scattering_texture, float4(i.uvw, 0));

				float3 deltaMultipleScatteringResult = texSample.rgb;
				float nu = texSample.a;

				return float4(mul(luminance_from_radiance, deltaMultipleScatteringResult) / RayleighPhaseFunction(nu), 0);
			}
			ENDCG
		}

	}
	Fallback Off
}
