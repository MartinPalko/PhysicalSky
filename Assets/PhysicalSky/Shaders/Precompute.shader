Shader "Hidden/Precompute"
{
	Properties {}
	SubShader
	{
		Pass
		{
			Cull Off ZWrite Off ZTest Always Fog{ Mode Off }

			CGPROGRAM
			#pragma target 5.0

			#pragma vertex vert_img
			#pragma fragment frag
			#pragma multi_compile TRANSMITTANCE SCATTERING IRRADIANCE
			
			#include "UnityCG.cginc"
			#include "PhysicalSkyCommon.cginc"

			RWTexture3D<float4> volumeTarget;
			float _volumeDepth;

			uniform IrradianceSpectrum _solar_irradiance;
			uniform Angle _sun_angular_radius;
			uniform Length _bottom_radius;
			uniform Length _top_radius;
			uniform Length _rayleigh_scale_height;
			uniform ScatteringSpectrum _rayleigh_scattering;
			uniform Length _mie_scale_height;
			uniform ScatteringSpectrum _mie_scattering;
			uniform ScatteringSpectrum _mie_extinction;
			uniform Number _mie_phase_function_g;
			uniform DimensionlessSpectrum _ground_albedo;
			uniform Number _mu_s_min;

			fixed4 frag (v2f_img i) : COLOR
			{
				fixed4 result = fixed4(1,1,1,1);

				AtmosphereParameters params;
				params.solar_irradiance = _solar_irradiance;
				params.sun_angular_radius = _sun_angular_radius;
				params.bottom_radius = _bottom_radius;
				params.top_radius = _top_radius;
				params.rayleigh_scale_height = _rayleigh_scale_height;
				params.rayleigh_scattering = _rayleigh_scattering;
				params.mie_scale_height = _mie_scale_height;
				params.mie_scattering = _mie_scattering;
				params.mie_extinction = _mie_extinction;
				params.mie_phase_function_g = _mie_phase_function_g;
				params.ground_albedo = _ground_albedo;
				params.mu_s_min = _mu_s_min;

#if TRANSMITTANCE
				result.rgb = ComputeTransmittanceToTopAtmosphereBoundaryTexture(params, i.uv);
#elif SCATTERING
				//int zz;
				//for (zz = 0; zz < _volumeDepth; zz++)
				//{
				//	float4 value = float4(i.uv.x, i.uv.y, zz / _volumeDepth, 1);
				//	volumeTarget[int3(i.pos.xy, zz)] = value;
				//}

				result.rgb = ComputeTransmittanceToTopAtmosphereBoundaryTexture(IN(AtmosphereParameters) atmosphere, IN(vec2) gl_frag_coord);
				discard;
#elif IRRADIANCE
				//result.rgb = ComputeTransmittanceToTopAtmosphereBoundaryTexture(params, i.uv);
#endif
				//return result.r == 0 ? float4(1,0,0,0) : float4(0,1,0,0);
				
				return result;

			}
			ENDCG
		}
	}
	Fallback Off
}
