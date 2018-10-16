#ifndef ATMOSPHERE_UNIFORMS
#define ATMOSPHERE_UNIFORMS

CBUFFER_START(atmosphere_parameters)
uniform float3 _sky_spectral_radiance_to_luminance;
uniform float3 _sun_spectral_radiance_to_luminance;

uniform IrradianceSpectrum _solar_irradiance;
uniform Angle _sun_angular_radius;
uniform Length _bottom_radius;
uniform Length _top_radius;
uniform float4 _rayleigh_density[3];
uniform ScatteringSpectrum _rayleigh_scattering;
uniform float4 _mie_density[3];
uniform ScatteringSpectrum _mie_scattering;
uniform ScatteringSpectrum _mie_extinction;
uniform Number _mie_phase_function_g;
uniform float4 _absorption_density[3];
uniform ScatteringSpectrum _absorption_extinction;
uniform DimensionlessSpectrum _ground_albedo;
uniform Number _mu_s_min;
CBUFFER_END

DensityProfile ArrayToDensityProfile(float4 a[3])
{
	DensityProfile profile;
	profile.layer0.width = a[0][0];
	profile.layer0.exp_term = a[0][1];
	profile.layer0.exp_scale = a[0][2];
	profile.layer0.linear_term = a[0][3];
	profile.layer0.constant_term = a[1][0];
	profile.layer1.width = a[1][1];
	profile.layer1.exp_term = a[1][2];
	profile.layer1.exp_scale = a[1][3];
	profile.layer1.linear_term = a[2][0];
	profile.layer1.constant_term = a[2][1];
	return profile;
}

AtmosphereParameters GetAtmosphereParameters()
{
	AtmosphereParameters params;

	params.sky_spectral_radiance_to_luminance = _sky_spectral_radiance_to_luminance;
	params.sun_spectral_radiance_to_luminance = _sun_spectral_radiance_to_luminance;

	params.solar_irradiance = _solar_irradiance;
	params.sun_angular_radius = _sun_angular_radius;
	params.bottom_radius = _bottom_radius;
	params.top_radius = _top_radius;
	params.rayleigh_density = ArrayToDensityProfile(_rayleigh_density);
	params.rayleigh_scattering = _rayleigh_scattering;
	params.mie_density = ArrayToDensityProfile(_mie_density);
	params.mie_scattering = _mie_scattering;
	params.mie_extinction = _mie_extinction;
	params.mie_phase_function_g = _mie_phase_function_g;
	params.absorption_density = ArrayToDensityProfile(_absorption_density);
	params.absorption_extinction = _absorption_extinction;
	params.ground_albedo = _ground_albedo;
	params.mu_s_min = _mu_s_min;
	return params;
}
#endif