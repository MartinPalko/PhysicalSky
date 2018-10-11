#ifndef ATMOSPHERE_UNIFORMS
#define ATMOSPHERE_UNIFORMS

uniform float3 _sky_spectral_radiance_to_luminance;
uniform float3 _sun_spectral_radiance_to_luminance;

uniform IrradianceSpectrum _solar_irradiance;
uniform Angle _sun_angular_radius;
uniform Length _bottom_radius;
uniform Length _top_radius;
uniform float _rayleigh_density0[5];
uniform float _rayleigh_density1[5];
uniform ScatteringSpectrum _rayleigh_scattering;
uniform float _mie_density0[5];
uniform float _mie_density1[5];
uniform ScatteringSpectrum _mie_scattering;
uniform ScatteringSpectrum _mie_extinction;
uniform Number _mie_phase_function_g;
uniform float _absorption_density0[5];
uniform float _absorption_density1[5];
uniform ScatteringSpectrum _absorption_extinction;
uniform DimensionlessSpectrum _ground_albedo;
uniform Number _mu_s_min;

DensityProfileLayer ArrayToDensityLayer(float a[5])
{
	DensityProfileLayer d;
	d.width = a[0];
	d.exp_term = a[1];
	d.exp_scale = a[2];
	d.linear_term = a[3];
	d.constant_term = a[4];
	return d;
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
	params.rayleigh_density.layer0 = ArrayToDensityLayer(_rayleigh_density0);
	params.rayleigh_density.layer1 = ArrayToDensityLayer(_rayleigh_density1);
	params.rayleigh_scattering = _rayleigh_scattering;
	params.mie_density.layer0 = ArrayToDensityLayer(_mie_density0);
	params.mie_density.layer1 = ArrayToDensityLayer(_mie_density1);
	params.mie_scattering = _mie_scattering;
	params.mie_extinction = _mie_extinction;
	params.mie_phase_function_g = _mie_phase_function_g;
	params.absorption_density.layer0 = ArrayToDensityLayer(_absorption_density0);
	params.absorption_density.layer1 = ArrayToDensityLayer(_absorption_density1);
	params.absorption_extinction = _absorption_extinction;
	params.ground_albedo = _ground_albedo;
	params.mu_s_min = _mu_s_min;
	return params;
}
#endif