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

AtmosphereParameters GetAtmosphereParameters()
{
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

	return params;
}