#ifndef PHYSICAL_SKY_CG
#define PHYSICAL_SKY_CG

#include "PhysicalSkyCommon.cginc"

uniform TransmittanceTexture transmittance_texture;
uniform ScatteringTexture scattering_texture;
uniform IrradianceTexture irradiance_texture;

float3 GetSolarLuminance()
{
	AtmosphereParameters params = GetAtmosphereParameters();
	return params.solar_irradiance /
		(PI * params.sun_angular_radius * params.sun_angular_radius) *
		params.sun_spectral_radiance_to_luminance;
}

float3 GetSkyLuminance(
	Position camera, Direction view_ray, Length shadow_length,
	Direction sun_direction, out DimensionlessSpectrum transmittance) 
{
	AtmosphereParameters params = GetAtmosphereParameters();
	return GetSkyRadiance(params, transmittance_texture,
		scattering_texture,
		camera, view_ray, shadow_length, sun_direction, transmittance) *
		params.sky_spectral_radiance_to_luminance;
}

float3 GetSkyLuminanceToPoint(
	Position camera, Position p, Length shadow_length,
	Direction sun_direction, out DimensionlessSpectrum transmittance)
{
	AtmosphereParameters params = GetAtmosphereParameters();
	return GetSkyRadianceToPoint(params, transmittance_texture,
		scattering_texture,
		camera, p, shadow_length, sun_direction, transmittance) *
		params.sky_spectral_radiance_to_luminance;
}

float3 GetSunAndSkyIlluminance(
	Position p, Direction normal, Direction sun_direction,
	out IrradianceSpectrum sky_irradiance) 
{
	AtmosphereParameters params = GetAtmosphereParameters();
	IrradianceSpectrum sun_irradiance = GetSunAndSkyIrradiance(
		params, transmittance_texture, irradiance_texture, p, normal,
		sun_direction, sky_irradiance);
	sky_irradiance *= params.sky_spectral_radiance_to_luminance;
	return sun_irradiance * params.sun_spectral_radiance_to_luminance;
}

#endif