/**
* Copyright (c) 2017 Eric Bruneton
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holders nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. in NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER in
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING in ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
* THE POSSIBILITY OF SUCH DAMAGE.
*
* Precomputed Atmospheric Scattering
* Copyright (c) 2008 INRIA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holders nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. in NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER in
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING in ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
* THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef PHYSICAL_SKY_COMMON
#define PHYSICAL_SKY_COMMON

#define NO_TEXTURE3D // Always enabled for now, since texture3d support in Unity is spotty.

#define TRANSMITTANCE_TEXTURE_WIDTH 256
#define TRANSMITTANCE_TEXTURE_HEIGHT 64
#define SCATTERING_TEXTURE_R_SIZE 32
#define SCATTERING_TEXTURE_MU_SIZE 128
#define SCATTERING_TEXTURE_MU_S_SIZE 32
#define SCATTERING_TEXTURE_NU_SIZE 8
#define SCATTERING_TEXTURE_WIDTH SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE
#ifdef NO_TEXTURE3D
	#define SCATTERING_TEXTURE_HEIGHT SCATTERING_TEXTURE_R_SIZE * SCATTERING_TEXTURE_MU_SIZE
	#define SCATTERING_TEXTURE_DEPTH 1
#else
	#define SCATTERING_TEXTURE_HEIGHT SCATTERING_TEXTURE_MU_SIZE
	#define SCATTERING_TEXTURE_DEPTH SCATTERING_TEXTURE_R_SIZE
#endif
#define IRRADIANCE_TEXTURE_WIDTH 64
#define IRRADIANCE_TEXTURE_HEIGHT 16

#define Length float
#define Wavelength float
#define Angle float
#define SolidAngle float
#define Power float
#define LuminousPower float

#define Number float
#define InverseLength float
#define Area float
#define Volume float
#define NumberDensity float
#define Irradiance float
#define Radiance float
#define SpectralPower float
#define SpectralIrradiance float
#define SpectralRadiance float
#define SpectralRadianceDensity float
#define ScatteringCoefficient float
#define InverseSolidAngle float
#define LuminousIntensity float
#define Luminance float
#define Illuminance float

// A generic function from Wavelength to some other type.
#define AbstractSpectrum float3
// A function from Wavelength to Number.
#define DimensionlessSpectrum float3
// A function from Wavelength to SpectralPower.
#define PowerSpectrum float3
// A function from Wavelength to SpectralIrradiance.
#define IrradianceSpectrum float3
// A function from Wavelength to SpectralRadiance.
#define RadianceSpectrum float3
// A function from Wavelength to SpectralRadianceDensity.
#define RadianceDensitySpectrum float3
// A function from Wavelength to ScaterringCoefficient.
#define ScatteringSpectrum float3

// A position in 3D (3 length values).
#define Position float3
// A unit direction vector in 3D (3 unitless values).
#define Direction float3
// A vector of 3 luminance values.
#define Luminance3 float3
// A vector of 3 illuminance values.
#define Illuminance3 float3

#define TransmittanceTexture sampler2D
#ifdef NO_TEXTURE3D
	#define AbstractScatteringTexture sampler2D
	#define ReducedScatteringTexture sampler2D
	#define ScatteringTexture sampler2D
	#define ScatteringDensityTexture sampler2D
#else
	#define AbstractScatteringTexture sampler3D
	#define ReducedScatteringTexture sampler3D
	#define ScatteringTexture sampler3D
	#define ScatteringDensityTexture sampler3D
#endif
#define IrradianceTexture sampler2D

static const Length m = 1.0; // Meter
static const Wavelength nm = 1.0; // NanoMeter
static const Angle rad = 1.0; // Radian
static const SolidAngle sr = 1.0;	// Steradian
static const Power watt = 1.0; // Watt
static const LuminousPower lm = 1.0; // Lumen

#ifndef PI
#define PI 3.1415926535897932386
#endif

static const Length km = 1000.0;
static const Area m2 = 1.0;
static const Volume m3 = 1.0;
static const Angle pi = 3.14159265358979323846;
static const Angle deg = 3.14159265358979323846 / 180.0;
static const Irradiance watt_per_square_meter = 1.0;
static const Radiance watt_per_square_meter_per_sr = 1.0;
static const SpectralIrradiance watt_per_square_meter_per_nm = 1.0;
static const SpectralRadiance watt_per_square_meter_per_sr_per_nm = 1.0;
static const SpectralRadianceDensity watt_per_cubic_meter_per_sr_per_nm = 1.0;
static const LuminousIntensity cd = 1.0;
static const LuminousIntensity kcd = 1000.0;
static const Luminance cd_per_square_meter = 1.0;
static const Luminance kcd_per_square_meter = 1000.0;

// An atmosphere layer of width 'width', and whose density is defined as
//   'exp_term' * exp('exp_scale' * h) + 'linear_term' * h + 'constant_term',
// clamped to [0,1], and where h is the altitude.
struct DensityProfileLayer
{
	Length width;
	Number exp_term;
	InverseLength exp_scale;
	InverseLength linear_term;
	Number constant_term;
};

// An atmosphere density profile made of several layers on top of each other
// (from bottom to top). The width of the last layer is ignored, i.e. it always
// extend to the top atmosphere boundary. The profile values vary between 0
// (null density) to 1 (maximum density).
struct DensityProfile 
{
	DensityProfileLayer layer0;
	DensityProfileLayer layer1;
};

//The atmosphere parameters are then defined by the following struct:
struct AtmosphereParameters 
{
	float3 sky_spectral_radiance_to_luminance;
	float3 sun_spectral_radiance_to_luminance;

	// The solar irradiance at the top of the atmosphere.
	IrradianceSpectrum solar_irradiance;
	// The sun's angular radius. Warning: the implementation uses approximations
	// that are valid only if this angle is smaller than 0.1 radians.
	Angle sun_angular_radius;
	// The distance between the planet center and the bottom of the atmosphere.
	Length bottom_radius;
	// The distance between the planet center and the top of the atmosphere.
	Length top_radius;
	// The density profile of air molecules, i.e. a function from altitude to
	// dimensionless values between 0 (null density) and 1 (maximum density).
	DensityProfile rayleigh_density;
	// The scattering coefficient of air molecules at the altitude where their
	// density is maximum (usually the bottom of the atmosphere), as a function of
	// wavelength. The scattering coefficient at altitude h is equal to
	// 'rayleigh_scattering' times 'rayleigh_density' at this altitude.
	ScatteringSpectrum rayleigh_scattering;
	// The density profile of aerosols, i.e. a function from altitude to
	// dimensionless values between 0 (null density) and 1 (maximum density).
	DensityProfile mie_density;
	// The scattering coefficient of aerosols at the altitude where their density
	// is maximum (usually the bottom of the atmosphere), as a function of
	// wavelength. The scattering coefficient at altitude h is equal to
	// 'mie_scattering' times 'mie_density' at this altitude.
	ScatteringSpectrum mie_scattering;
	// The extinction coefficient of aerosols at the altitude where their density
	// is maximum (usually the bottom of the atmosphere), as a function of
	// wavelength. The extinction coefficient at altitude h is equal to
	// 'mie_extinction' times 'mie_density' at this altitude.
	ScatteringSpectrum mie_extinction;
	// The asymetry parameter for the Cornette-Shanks phase function for the
	// aerosols.
	Number mie_phase_function_g;
	// The density profile of air molecules that absorb light (e.g. ozone), i.e.
	// a function from altitude to dimensionless values between 0 (null density)
	// and 1 (maximum density).
	DensityProfile absorption_density;
	// The extinction coefficient of molecules that absorb light (e.g. ozone) at
	// the altitude where their density is maximum, as a function of wavelength.
	// The extinction coefficient at altitude h is equal to
	// 'absorption_extinction' times 'absorption_density' at this altitude.
	ScatteringSpectrum absorption_extinction;
	// The average albedo of the ground.
	DimensionlessSpectrum ground_albedo;
	// The cosine of the maximum Sun zenith angle for which atmospheric scattering
	// must be precomputed (for maximum precision, use the smallest Sun zenith
	// angle yielding negligible sky light radiance values. For instance, for the
	// Earth case, 102 degrees is a good choice - yielding mu_s_min = -0.2).
	Number mu_s_min;
};

#include "AtmosphereUniforms.cginc"

float4 scatteringTextureSample(ScatteringTexture s, float3 uvw)
{
#ifdef NO_TEXTURE3D
	float tex_coord_z = uvw.z * Number(SCATTERING_TEXTURE_R_SIZE - 1);
	float tex_z = floor(tex_coord_z);
	float l = tex_coord_z - tex_z;
	float2 uv0 = float2(uvw.x, (tex_z + uvw.y) / Number(SCATTERING_TEXTURE_R_SIZE));
	float2 uv1 = float2(uvw.x, (tex_z + 1.0 + uvw.y) / Number(SCATTERING_TEXTURE_R_SIZE));
	return lerp(tex2Dlod(s, float4(uv0, 0, 0)), tex2Dlod(s, float4(uv1, 0, 0)), l);
#else
	return tex3Dlod(s, float4(uvw, 0));
#endif
}

float4 scatteringTextureSample(ScatteringTexture s, float4 uvwz)
{
	float tex_coord_x = uvwz.x * Number(SCATTERING_TEXTURE_NU_SIZE - 1);
	float tex_x = floor(tex_coord_x);
	float l = tex_coord_x - tex_x;
	float3 uvw0 = float3((tex_x + uvwz.y) / Number(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
	float3 uvw1 = float3((tex_x + 1.0 + uvwz.y) / Number(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
	return lerp(scatteringTextureSample(s, uvw0), scatteringTextureSample(s, uvw1), l);
}

Number ClampCosine(Number mu)
{
	return clamp(mu, Number(-1.0), Number(1.0));
}

Length ClampDistance(Length d)
{
	return max(d, 0.0 * m);
}

Length ClampRadius(in AtmosphereParameters atmosphere, Length r) 
{
	return clamp(r, atmosphere.bottom_radius, atmosphere.top_radius);
}

Length SafeSqrt(Area a) 
{
	return sqrt(max(a, 0.0 * m2));
}

RadianceSpectrum GetSolarRadiance(AtmosphereParameters atmosphere) 
{
	return atmosphere.solar_irradiance /
		(PI * atmosphere.sun_angular_radius * atmosphere.sun_angular_radius);
}

Length DistanceToTopAtmosphereBoundary(in AtmosphereParameters atmosphere,	Length r, Number mu) 
{
	//assert(r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	Area discriminant = r * r * (mu * mu - 1.0) +
		atmosphere.top_radius * atmosphere.top_radius;
	return ClampDistance(-r * mu + SafeSqrt(discriminant));
}

Length DistanceToBottomAtmosphereBoundary(in AtmosphereParameters atmosphere, Length r, Number mu)
{
	//assert(r >= atmosphere.bottom_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	Area discriminant = r * r * (mu * mu - 1.0) +
		atmosphere.bottom_radius * atmosphere.bottom_radius;
	return ClampDistance(-r * mu - SafeSqrt(discriminant));
}

bool RayIntersectsGround(in AtmosphereParameters atmosphere, Length r, Number mu) 
{
	//assert(r >= atmosphere.bottom_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	return mu < 0.0 && r * r * (mu * mu - 1.0) +
		atmosphere.bottom_radius * atmosphere.bottom_radius >= 0.0 * m2;
}

Number GetLayerDensity(in DensityProfileLayer layer, Length altitude) 
{
	Number density = layer.exp_term * exp(layer.exp_scale * altitude) +
		layer.linear_term * altitude + layer.constant_term;
	return clamp(density, Number(0.0), Number(1.0));
}

Number GetProfileDensity(in DensityProfile profile, Length altitude) 
{
	return altitude < profile.layer0.width ?
		GetLayerDensity(profile.layer0, altitude) :
		GetLayerDensity(profile.layer1, altitude);
}

Length ComputeOpticalLengthToTopAtmosphereBoundary(in AtmosphereParameters atmosphere, in DensityProfile profile, Length r, Number mu)
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	// Number of intervals for the numerical integration.
	const int SAMPLE_COUNT = 500;
	// The integration step, i.e. the length of each integration interval.
	Length dx =	DistanceToTopAtmosphereBoundary(atmosphere, r, mu) / Number(SAMPLE_COUNT);
	// Integration loop.
	Length result = 0.0 * m;
	for (int i = 0; i <= SAMPLE_COUNT; ++i)
	{
		Length d_i = Number(i) * dx;
		// Distance between the current sample _point and the planet center.
		Length r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
		// Number density at the current sample _point (divided by the number density
		// at the bottom of the atmosphere, yielding a dimensionless number).
		Number y_i = GetProfileDensity(profile, r_i - atmosphere.bottom_radius);
		// Sample weight (from the trapezoidal rule).
		Number weight_i = i == 0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
		result += y_i * weight_i * dx;
	}
	return result;
}

DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundary(in AtmosphereParameters atmosphere, Length r, Number mu)
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	return exp(-(
		atmosphere.rayleigh_scattering *
		ComputeOpticalLengthToTopAtmosphereBoundary(
			atmosphere, atmosphere.rayleigh_density, r, mu) +
		atmosphere.mie_extinction *
		ComputeOpticalLengthToTopAtmosphereBoundary(
			atmosphere, atmosphere.mie_density, r, mu) +
      atmosphere.absorption_extinction *
          ComputeOpticalLengthToTopAtmosphereBoundary(
              atmosphere, atmosphere.absorption_density, r, mu)));
}

Number GetTextureCoordFromUnitRange(Number x, int texture_size)
{
	return 0.5 / Number(texture_size) + x * (1.0 - 1.0 / Number(texture_size));
}

Number GetUnitRangeFromTextureCoord(Number u, int texture_size)
{
	return (u - 0.5 / Number(texture_size)) / (1.0 - 1.0 / Number(texture_size));
}

float2 GetTransmittanceTextureUvFromRMu(in AtmosphereParameters atmosphere, Length r, Number mu) 
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	// Distance to top atmosphere boundary for a horizontal ray at ground level.
	Length H = sqrt(atmosphere.top_radius * atmosphere.top_radius -	atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the horizon.
	Length rho = SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
	// and maximum values over all mu - obtained for (r,1) and (r,mu_horizon).
	Length d = DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
	Length d_min = atmosphere.top_radius - r;
	Length d_max = rho + H;
	Number x_mu = (d - d_min) / (d_max - d_min);
	Number x_r = rho / H;
	return float2(GetTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_TEXTURE_WIDTH), GetTextureCoordFromUnitRange(x_r, TRANSMITTANCE_TEXTURE_HEIGHT));
}

void GetRMuFromTransmittanceTextureUv(in AtmosphereParameters atmosphere, in float2 uv, out Length r, out Number mu) 
{
	//assert(uv.x >= 0.0 && uv.x <= 1.0);
	//assert(uv.y >= 0.0 && uv.y <= 1.0);
	Number x_mu = GetUnitRangeFromTextureCoord(uv.x, TRANSMITTANCE_TEXTURE_WIDTH);
	Number x_r = GetUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT);
	// Distance to top atmosphere boundary for a horizontal ray at ground level.
	Length H = sqrt(atmosphere.top_radius * atmosphere.top_radius -	atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the horizon, from which we can compute r:
	Length rho = H * x_r;
	r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
	// and maximum values over all mu - obtained for (r,1) and (r,mu_horizon) -
	// from which we can recover mu:
	Length d_min = atmosphere.top_radius - r;
	Length d_max = rho + H;
	Length d = d_min + x_mu * (d_max - d_min);
	mu = d == 0.0 * m ? Number(1.0) : (H * H - rho * rho - d * d) / (2.0 * r * d);
	mu = ClampCosine(mu);
}

DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundaryTexture(in AtmosphereParameters atmosphere, in float2 frag_coord)
{
	// Don't have to scale by texture size, since we now just pass in normalized values.
	const float2 TRANSMITTANCE_TEXTURE_SIZE =	float2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
	Length r;
	Number mu;
	GetRMuFromTransmittanceTextureUv(atmosphere, frag_coord / TRANSMITTANCE_TEXTURE_SIZE, r, mu);
	return ComputeTransmittanceToTopAtmosphereBoundary(atmosphere, r, mu);
}

DimensionlessSpectrum GetTransmittanceToTopAtmosphereBoundary(in AtmosphereParameters atmosphere, in TransmittanceTexture transmittance_texture, Length r, Number mu)
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	float2 uv = GetTransmittanceTextureUvFromRMu(atmosphere, r, mu);
	return DimensionlessSpectrum(tex2Dlod(transmittance_texture, float4(uv, 0, 0)).rgb);
}

DimensionlessSpectrum GetTransmittance(in AtmosphereParameters atmosphere, in TransmittanceTexture transmittance_texture, Length r, Number mu, Length d, bool ray_r_mu_intersects_ground) 
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	//assert(d >= 0.0 * m);

	Length r_d = ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
	Number mu_d = ClampCosine((r * mu + d) / r_d);

	if (ray_r_mu_intersects_ground)
	{
		return min(
			GetTransmittanceToTopAtmosphereBoundary(
				atmosphere, transmittance_texture, r_d, -mu_d) /
			GetTransmittanceToTopAtmosphereBoundary(
				atmosphere, transmittance_texture, r, -mu),
			DimensionlessSpectrum(1.0, 1.0, 1.0));
	}
	else 
	{
		return min(
			GetTransmittanceToTopAtmosphereBoundary(
				atmosphere, transmittance_texture, r, mu) /
			GetTransmittanceToTopAtmosphereBoundary(
				atmosphere, transmittance_texture, r_d, mu_d),
			DimensionlessSpectrum(1.0, 1.0, 1.0));
	}
}

DimensionlessSpectrum GetTransmittanceToSun(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	Length r, Number mu_s)
{
	Number sin_theta_h = atmosphere.bottom_radius / r;
	Number cos_theta_h = -sqrt(max(1.0 - sin_theta_h * sin_theta_h, 0.0));
	return GetTransmittanceToTopAtmosphereBoundary(
		atmosphere, transmittance_texture, r, mu_s) *
		smoothstep(-sin_theta_h * atmosphere.sun_angular_radius / rad,
			sin_theta_h * atmosphere.sun_angular_radius / rad,
			mu_s - cos_theta_h);
}

void ComputeSingleScatteringIntegrand(in AtmosphereParameters atmosphere, in TransmittanceTexture transmittance_texture, Length r, Number mu, Number mu_s, Number nu, Length d, bool ray_r_mu_intersects_ground, out DimensionlessSpectrum rayleigh, out DimensionlessSpectrum mie)
{
	Length r_d = ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
	Number mu_s_d = ClampCosine((r * mu_s + d * nu) / r_d);

	DimensionlessSpectrum transmittance =
		GetTransmittance(
			atmosphere, transmittance_texture, r, mu, d,
			ray_r_mu_intersects_ground) *
		GetTransmittanceToSun(
			atmosphere, transmittance_texture, r_d, mu_s_d);
	rayleigh = transmittance * GetProfileDensity(
		atmosphere.rayleigh_density, r_d - atmosphere.bottom_radius);
	mie = transmittance * GetProfileDensity(
		atmosphere.mie_density, r_d - atmosphere.bottom_radius);
}

Length DistanceToNearestAtmosphereBoundary(in AtmosphereParameters atmosphere,	Length r, Number mu, bool ray_r_mu_intersects_ground) 
{
	if (ray_r_mu_intersects_ground) 
	{
		return DistanceToBottomAtmosphereBoundary(atmosphere, r, mu);
	}
	else
	{
		return DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
	}
}

void ComputeSingleScattering(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground,
	out IrradianceSpectrum rayleigh, out IrradianceSpectrum mie) 
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	//assert(mu_s >= -1.0 && mu_s <= 1.0);
	//assert(nu >= -1.0 && nu <= 1.0);

	// Number of intervals for the numerical integration.
	const int SAMPLE_COUNT = 50;
	// The integration step, i.e. the length of each integration interval.
	Length dx =	DistanceToNearestAtmosphereBoundary(atmosphere, r, mu, ray_r_mu_intersects_ground) / Number(SAMPLE_COUNT);
	// Integration loop.
	DimensionlessSpectrum rayleigh_sum = DimensionlessSpectrum(0.0, 0.0, 0.0);
	DimensionlessSpectrum mie_sum = DimensionlessSpectrum(0.0, 0.0, 0.0);
	for (int i = 0; i <= SAMPLE_COUNT; ++i)
	{
		Length d_i = Number(i) * dx;
		// The Rayleigh and Mie single scattering at the current sample _point.
		DimensionlessSpectrum rayleigh_i;
		DimensionlessSpectrum mie_i;
		ComputeSingleScatteringIntegrand(atmosphere, transmittance_texture,
			r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground, rayleigh_i, mie_i);
		// Sample weight (from the trapezoidal rule).
		Number weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
		rayleigh_sum += rayleigh_i * weight_i;
		mie_sum += mie_i * weight_i;
	}
	rayleigh = rayleigh_sum * dx * atmosphere.solar_irradiance * atmosphere.rayleigh_scattering;
	mie = mie_sum * dx * atmosphere.solar_irradiance * atmosphere.mie_scattering;
}

InverseSolidAngle RayleighPhaseFunction(Number nu)
{
	InverseSolidAngle k = 3.0 / (16.0 * PI * sr);
	return k * (1.0 + nu * nu);
}

InverseSolidAngle MiePhaseFunction(Number g, Number nu) 
{
	InverseSolidAngle k = 3.0 / (8.0 * PI * sr) * (1.0 - g * g) / (2.0 + g * g);
	return k * (1.0 + nu * nu) / pow(abs(1.0 + g * g - 2.0 * g * nu), 1.5);
}

float4 GetScatteringTextureUvwzFromRMuMuSNu(in AtmosphereParameters atmosphere, Length r, Number mu, Number mu_s, Number nu, bool ray_r_mu_intersects_ground) 
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	//assert(mu_s >= -1.0 && mu_s <= 1.0);
	//assert(nu >= -1.0 && nu <= 1.0);

	// Distance to top atmosphere boundary for a horizontal ray at ground level.
	Length H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
		atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the horizon.
	Length rho =
		SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);
	Number u_r = GetTextureCoordFromUnitRange(rho / H, SCATTERING_TEXTURE_R_SIZE);

	// Discriminant of the quadratic equation for the intersections of the ray
	// (r,mu) with the ground (see RayIntersectsGround).
	Length r_mu = r * mu;
	Area discriminant =
		r_mu * r_mu - r * r + atmosphere.bottom_radius * atmosphere.bottom_radius;
	Number u_mu;
	if (ray_r_mu_intersects_ground) 
	{
		// Distance to the ground for the ray (r,mu), and its minimum and maximum
		// values over all mu - obtained for (r,-1) and (r,mu_horizon).
		Length d = -r_mu - SafeSqrt(discriminant);
		Length d_min = r - atmosphere.bottom_radius;
		Length d_max = rho;
		u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(d_max == d_min ? 0.0 :
			(d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
	}
	else
	{
		// Distance to the top atmosphere boundary for the ray (r,mu), and its
		// minimum and maximum values over all mu - obtained for (r,1) and
		// (r,mu_horizon).
		Length d = -r_mu + SafeSqrt(discriminant + H * H);
		Length d_min = atmosphere.top_radius - r;
		Length d_max = rho + H;
		u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
			(d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
	}

	Length d = DistanceToTopAtmosphereBoundary(atmosphere, atmosphere.bottom_radius, mu_s);
	Length d_min = atmosphere.top_radius - atmosphere.bottom_radius;
	Length d_max = H;
	Number a = (d - d_min) / (d_max - d_min);
	Number A = -2.0 * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
	Number u_mu_s = GetTextureCoordFromUnitRange(max(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_S_SIZE);
	Number u_nu = (nu + 1.0) / 2.0;
	return float4(u_nu, u_mu_s, u_mu, u_r);
}

void GetRMuMuSNuFromScatteringTextureUvwz(in AtmosphereParameters atmosphere,
	in float4 uvwz, out Length r, out Number mu, out Number mu_s,
	out Number nu, out bool ray_r_mu_intersects_ground)
{
	// Distance to top atmosphere boundary for a horizontal ray at ground level.
	Length H = sqrt(atmosphere.top_radius * atmosphere.top_radius -	atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the horizon.
	Length rho =
		H * GetUnitRangeFromTextureCoord(uvwz.w, SCATTERING_TEXTURE_R_SIZE);
	r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);

	if (uvwz.z < 0.5)
	{
		// Distance to the ground for the ray (r,mu), and its minimum and maximum
		// values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
		// we can recover mu:
		Length d_min = r - atmosphere.bottom_radius;
		Length d_max = rho;
		Length d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
			1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE / 2);
		mu = d == 0.0 * m ? Number(-1.0) :
			ClampCosine(-(rho * rho + d * d) / (2.0 * r * d));
		ray_r_mu_intersects_ground = true;
	}
	else 
	{
		// Distance to the top atmosphere boundary for the ray (r,mu), and its
		// minimum and maximum values over all mu - obtained for (r,1) and
		// (r,mu_horizon) - from which we can recover mu:
		Length d_min = atmosphere.top_radius - r;
		Length d_max = rho + H;
		Length d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
			2.0 * uvwz.z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2);
		mu = d == 0.0 * m ? Number(1.0) :
			ClampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
		ray_r_mu_intersects_ground = false;
	}

	Number x_mu_s =
		GetUnitRangeFromTextureCoord(uvwz.y, SCATTERING_TEXTURE_MU_S_SIZE);
	Length d_min = atmosphere.top_radius - atmosphere.bottom_radius;
	Length d_max = H;
	Number A =
		-2.0 * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
	Number a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
	Length d = d_min + min(a, A) * (d_max - d_min);
	mu_s = d == 0.0 * m ? Number(1.0) :
		ClampCosine((H * H - d * d) / (2.0 * atmosphere.bottom_radius * d));

	nu = ClampCosine(uvwz.x * 2.0 - 1.0);
}

void GetRMuMuSNuFromScatteringTextureFragCoord(
	in AtmosphereParameters atmosphere, in float3 frag_coord,
	out Length r, out Number mu, out Number mu_s, out Number nu,
	out bool ray_r_mu_intersects_ground)
{
	const float4 SCATTERING_TEXTURE_SIZE = float4(
		SCATTERING_TEXTURE_NU_SIZE - 1,
		SCATTERING_TEXTURE_MU_S_SIZE,
		SCATTERING_TEXTURE_MU_SIZE,
#ifdef NO_TEXTURE3D
		SCATTERING_TEXTURE_R_SIZE - 1
#else
		SCATTERING_TEXTURE_R_SIZE
#endif
		);
	
	float4 uvwz;
	uvwz.x = floor(frag_coord.x / Number(SCATTERING_TEXTURE_MU_S_SIZE));
	uvwz.y = fmod(frag_coord.x, Number(SCATTERING_TEXTURE_MU_S_SIZE));
	
#ifdef NO_TEXTURE3D
	uvwz.w = floor(frag_coord.y / Number(SCATTERING_TEXTURE_MU_SIZE));
	uvwz.z = fmod(frag_coord.y, Number(SCATTERING_TEXTURE_MU_SIZE));
#else
	uvwz.z = frag_coord.y;
	uvwz.w = frag_coord.z;
#endif

	uvwz /= SCATTERING_TEXTURE_SIZE;

	GetRMuMuSNuFromScatteringTextureUvwz(
		atmosphere, uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	// Clamp nu to its valid range of values, given mu and mu_s.
	nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
		mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}

void ComputeSingleScatteringTexture(in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture, in float3 frag_coord,
	out IrradianceSpectrum rayleigh, out IrradianceSpectrum mie)
{
	Length r;
	Number mu;
	Number mu_s;
	Number nu;
	bool ray_r_mu_intersects_ground;
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
		r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	ComputeSingleScattering(atmosphere, transmittance_texture,
		r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
}

AbstractSpectrum GetScattering(
	in AtmosphereParameters atmosphere,
	in AbstractScatteringTexture scattering_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground)
{
	float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(
		atmosphere, r, mu, mu_s, nu, ray_r_mu_intersects_ground);

	return AbstractSpectrum(scatteringTextureSample(scattering_texture, uvwz).rgb);
}

RadianceSpectrum GetScattering(
	in AtmosphereParameters atmosphere,
	in ReducedScatteringTexture single_rayleigh_scattering_texture,
	in ReducedScatteringTexture single_mie_scattering_texture,
	in ScatteringTexture multiple_scattering_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground,
	int scattering_order)
{
	if (scattering_order == 1) 
	{
		IrradianceSpectrum rayleigh = GetScattering(
			atmosphere, single_rayleigh_scattering_texture, r, mu, mu_s, nu,
			ray_r_mu_intersects_ground);
		IrradianceSpectrum mie = GetScattering(
			atmosphere, single_mie_scattering_texture, r, mu, mu_s, nu,
			ray_r_mu_intersects_ground);
		return rayleigh * RayleighPhaseFunction(nu) +
			mie * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
	}
	else
	{
		return GetScattering(
			atmosphere, multiple_scattering_texture, r, mu, mu_s, nu,
			ray_r_mu_intersects_ground);
	}
}

IrradianceSpectrum GetIrradiance(
	in AtmosphereParameters atmosphere,
	in IrradianceTexture irradiance_texture,
	Length r, Number mu_s);


RadianceDensitySpectrum ComputeScatteringDensity(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	in ReducedScatteringTexture single_rayleigh_scattering_texture,
	in ReducedScatteringTexture single_mie_scattering_texture,
	in ScatteringTexture multiple_scattering_texture,
	in IrradianceTexture irradiance_texture,
	Length r, Number mu, Number mu_s, Number nu, int scattering_order)
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	//assert(mu_s >= -1.0 && mu_s <= 1.0);
	//assert(nu >= -1.0 && nu <= 1.0);
	//assert(scattering_order >= 2);

	// Compute unit direction vectors for the zenith, the view direction omega and
	// and the sun direction omega_s, such that the cosine of the view-zenith
	// angle is mu, the cosine of the sun-zenith angle is mu_s, and the cosine of
	// the view-sun angle is nu. The goal is to simplify computations below.
	float3 zenith_direction = float3(0.0, 0.0, 1.0);
	float3 omega = float3(sqrt(1.0 - mu * mu), 0.0, mu);
	Number sun_dir_x = omega.x == 0.0 ? 0.0 : (nu - mu * mu_s) / omega.x;
	Number sun_dir_y = sqrt(max(1.0 - sun_dir_x * sun_dir_x - mu_s * mu_s, 0.0));
	float3 omega_s = float3(sun_dir_x, sun_dir_y, mu_s);

	const int SAMPLE_COUNT = 16;
	const Angle dphi = pi / Number(SAMPLE_COUNT);
	const Angle dtheta = pi / Number(SAMPLE_COUNT);
	RadianceDensitySpectrum rayleigh_mie =
		RadianceDensitySpectrum(0.0, 0.0, 0.0);

	// Nested loops for the integral over all the incident directions omega_i.
	for (int l = 0; l < SAMPLE_COUNT; ++l)
	{
		Angle theta = (Number(l) + 0.5) * dtheta;
		Number cos_theta = cos(theta);
		Number sin_theta = sin(theta);
		bool ray_r_theta_intersects_ground =
			RayIntersectsGround(atmosphere, r, cos_theta);

		// The distance and transmittance to the ground only depend on theta, so we
		// can compute them in the outer loop for efficiency.
		Length distance_to_ground = 0.0 * m;
		DimensionlessSpectrum transmittance_to_ground = DimensionlessSpectrum(0.0, 0.0, 0.0);
		DimensionlessSpectrum ground_albedo = DimensionlessSpectrum(0.0, 0.0, 0.0);
		if (ray_r_theta_intersects_ground)
		{
			distance_to_ground =
				DistanceToBottomAtmosphereBoundary(atmosphere, r, cos_theta);
			transmittance_to_ground =
				GetTransmittance(atmosphere, transmittance_texture, r, cos_theta,
					distance_to_ground, true /* ray_intersects_ground */);
			ground_albedo = atmosphere.ground_albedo;
		}

		for (int m = 0; m < 2 * SAMPLE_COUNT; ++m)
		{
			Angle phi = (Number(m) + 0.5) * dphi;
			float3 omega_i =
				float3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
			SolidAngle domega_i = (dtheta / rad) * (dphi / rad) * sin(theta) * sr;

			// The radiance L_i arriving from direction omega_i after n-1 bounces is
			// the sum of a term given by the precomputed scattering texture for the
			// (n-1)-th order:
			Number nu1 = dot(omega_s, omega_i);
			RadianceSpectrum incident_radiance = GetScattering(atmosphere,
				single_rayleigh_scattering_texture, single_mie_scattering_texture,
				multiple_scattering_texture, r, omega_i.z, mu_s, nu1,
				ray_r_theta_intersects_ground, scattering_order - 1);

			// and of the contribution from the light paths with n-1 bounces and whose
			// last bounce is on the ground. This contribution is the product of the
			// transmittance to the ground, the ground albedo, the ground BRDF, and
			// the irradiance received on the ground after n-2 bounces.
			float3 ground_normal =
				normalize(zenith_direction * r + omega_i * distance_to_ground);
			IrradianceSpectrum ground_irradiance = GetIrradiance(
				atmosphere, irradiance_texture, atmosphere.bottom_radius,
				dot(ground_normal, omega_s));
			incident_radiance += transmittance_to_ground *
				ground_albedo * (1.0 / (PI * sr)) * ground_irradiance;

			// The radiance finally scattered from direction omega_i towards direction
			// -omega is the product of the incident radiance, the scattering
			// coefficient, and the phase function for directions omega and omega_i
			// (all this summed over all particle types, i.e. Rayleigh and Mie).
			Number nu2 = dot(omega, omega_i);
			Number rayleigh_density = GetProfileDensity(
				atmosphere.rayleigh_density, r - atmosphere.bottom_radius);
			Number mie_density = GetProfileDensity(
				atmosphere.mie_density, r - atmosphere.bottom_radius);
			rayleigh_mie += incident_radiance * (
				atmosphere.rayleigh_scattering * rayleigh_density *
				RayleighPhaseFunction(nu2) +
				atmosphere.mie_scattering * mie_density *
				MiePhaseFunction(atmosphere.mie_phase_function_g, nu2)) *
				domega_i;
		}
	}
	return rayleigh_mie;
}

RadianceSpectrum ComputeMultipleScattering(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	in ScatteringDensityTexture scattering_density_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground)
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu >= -1.0 && mu <= 1.0);
	//assert(mu_s >= -1.0 && mu_s <= 1.0);
	//assert(nu >= -1.0 && nu <= 1.0);

	// Number of intervals for the numerical integration.
	const int SAMPLE_COUNT = 50;
	// The integration step, i.e. the length of each integration interval.
	Length dx =
		DistanceToNearestAtmosphereBoundary(
			atmosphere, r, mu, ray_r_mu_intersects_ground) /
		Number(SAMPLE_COUNT);
	// Integration loop.
	RadianceSpectrum rayleigh_mie_sum =
		RadianceSpectrum(0.0, 0.0, 0.0);
	for (int i = 0; i <= SAMPLE_COUNT; ++i)
	{
		Length d_i = Number(i) * dx;

		// The r, mu and mu_s parameters at the current integration _point (see the
		// single scattering section for a detailed explanation).
		Length r_i =
			ClampRadius(atmosphere, sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r));
		Number mu_i = ClampCosine((r * mu + d_i) / r_i);
		Number mu_s_i = ClampCosine((r * mu_s + d_i * nu) / r_i);

		// The Rayleigh and Mie multiple scattering at the current sample _point.
		RadianceSpectrum rayleigh_mie_i =
			GetScattering(
				atmosphere, scattering_density_texture, r_i, mu_i, mu_s_i, nu,
				ray_r_mu_intersects_ground) *
			GetTransmittance(
				atmosphere, transmittance_texture, r, mu, d_i,
				ray_r_mu_intersects_ground) *
			dx;
		// Sample weight (from the trapezoidal rule).
		Number weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
		rayleigh_mie_sum += rayleigh_mie_i * weight_i;
	}
	return rayleigh_mie_sum;
}

RadianceDensitySpectrum ComputeScatteringDensityTexture(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	in ReducedScatteringTexture single_rayleigh_scattering_texture,
	in ReducedScatteringTexture single_mie_scattering_texture,
	in ScatteringTexture multiple_scattering_texture,
	in IrradianceTexture irradiance_texture,
	in float3 frag_coord, int scattering_order) 
{
	Length r;
	Number mu;
	Number mu_s;
	Number nu;
	bool ray_r_mu_intersects_ground;
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
		r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	return ComputeScatteringDensity(atmosphere, transmittance_texture,
		single_rayleigh_scattering_texture, single_mie_scattering_texture,
		multiple_scattering_texture, irradiance_texture, r, mu, mu_s, nu,
		scattering_order);
}

RadianceSpectrum ComputeMultipleScatteringTexture(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	in ScatteringDensityTexture scattering_density_texture,
	in float3 frag_coord, out Number nu)
{
	Length r;
	Number mu;
	Number mu_s;
	bool ray_r_mu_intersects_ground;
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
		r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	return ComputeMultipleScattering(atmosphere, transmittance_texture,
		scattering_density_texture, r, mu, mu_s, nu,
		ray_r_mu_intersects_ground);
}

IrradianceSpectrum ComputeDirectIrradiance(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	Length r, Number mu_s) 
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu_s >= -1.0 && mu_s <= 1.0);

	Number alpha_s = atmosphere.sun_angular_radius / rad;
	// Approximate average of the cosine factor mu_s over the visible fraction of
	// the Sun disc.
	Number average_cosine_factor =
		mu_s < -alpha_s ? 0.0 : (mu_s > alpha_s ? mu_s :
		(mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s));

	return atmosphere.solar_irradiance *
		GetTransmittanceToTopAtmosphereBoundary(
			atmosphere, transmittance_texture, r, mu_s) * average_cosine_factor;
}

IrradianceSpectrum ComputeIndirectIrradiance(
	in AtmosphereParameters atmosphere,
	in ReducedScatteringTexture single_rayleigh_scattering_texture,
	in ReducedScatteringTexture single_mie_scattering_texture,
	in ScatteringTexture multiple_scattering_texture,
	Length r, Number mu_s, int scattering_order) 
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu_s >= -1.0 && mu_s <= 1.0);
	//assert(scattering_order >= 1);

	const int SAMPLE_COUNT = 32;
	const Angle dphi = pi / Number(SAMPLE_COUNT);
	const Angle dtheta = pi / Number(SAMPLE_COUNT);

	IrradianceSpectrum result =
		IrradianceSpectrum(0.0, 0.0, 0.0);
	float3 omega_s = float3(sqrt(1.0 - mu_s * mu_s), 0.0, mu_s);
	for (int j = 0; j < SAMPLE_COUNT / 2; ++j)
	{
		Angle theta = (Number(j) + 0.5) * dtheta;
		for (int i = 0; i < 2 * SAMPLE_COUNT; ++i) 
		{
			Angle phi = (Number(i) + 0.5) * dphi;
			float3 omega =
				float3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
			SolidAngle domega = (dtheta / rad) * (dphi / rad) * sin(theta) * sr;

			Number nu = dot(omega, omega_s);
			result += GetScattering(atmosphere, single_rayleigh_scattering_texture,
				single_mie_scattering_texture, multiple_scattering_texture,
				r, omega.z, mu_s, nu, false /* ray_r_theta_intersects_ground */,
				scattering_order) *
				omega.z * domega;
		}
	}
	return result;
}

float2 GetIrradianceTextureUvFromRMuS(in AtmosphereParameters atmosphere, Length r, Number mu_s) 
{
	//assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	//assert(mu_s >= -1.0 && mu_s <= 1.0);
	Number x_r = (r - atmosphere.bottom_radius) / (atmosphere.top_radius - atmosphere.bottom_radius);
	Number x_mu_s = mu_s * 0.5 + 0.5;
	return float2(GetTextureCoordFromUnitRange(x_mu_s, IRRADIANCE_TEXTURE_WIDTH),	GetTextureCoordFromUnitRange(x_r, IRRADIANCE_TEXTURE_HEIGHT));
}

void GetRMuSFromIrradianceTextureUv(in AtmosphereParameters atmosphere, in float2 uv, out Length r, out Number mu_s)
{
	//assert(uv.x >= 0.0 && uv.x <= 1.0);
	//assert(uv.y >= 0.0 && uv.y <= 1.0);
	Number x_mu_s = GetUnitRangeFromTextureCoord(uv.x, IRRADIANCE_TEXTURE_WIDTH);
	Number x_r = GetUnitRangeFromTextureCoord(uv.y, IRRADIANCE_TEXTURE_HEIGHT);
	r = atmosphere.bottom_radius + x_r * (atmosphere.top_radius - atmosphere.bottom_radius);
	mu_s = ClampCosine(2.0 * x_mu_s - 1.0);
}

#define IRRADIANCE_TEXTURE_SIZE float2(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT)

IrradianceSpectrum ComputeDirectIrradianceTexture(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	in float2 frag_coord) 
{
	Length r;
	Number mu_s;
	GetRMuSFromIrradianceTextureUv(atmosphere, frag_coord / IRRADIANCE_TEXTURE_SIZE, r, mu_s);
	return ComputeDirectIrradiance(atmosphere, transmittance_texture, r, mu_s);
}

IrradianceSpectrum ComputeIndirectIrradianceTexture(
	in AtmosphereParameters atmosphere,
	in ReducedScatteringTexture single_rayleigh_scattering_texture,
	in ReducedScatteringTexture single_mie_scattering_texture,
	in ScatteringTexture multiple_scattering_texture,
	in float2 frag_coord, int scattering_order)
{
	Length r;
	Number mu_s;
	GetRMuSFromIrradianceTextureUv(
		atmosphere, frag_coord / IRRADIANCE_TEXTURE_SIZE, r, mu_s);
	return ComputeIndirectIrradiance(atmosphere, single_rayleigh_scattering_texture, single_mie_scattering_texture,	multiple_scattering_texture, r, mu_s, scattering_order);
}

IrradianceSpectrum GetIrradiance(
	in AtmosphereParameters atmosphere,
	in IrradianceTexture irradiance_texture,
	Length r, Number mu_s) 
{
	float2 uv = GetIrradianceTextureUvFromRMuS(atmosphere, r, mu_s);
	return IrradianceSpectrum(tex2Dlod(irradiance_texture, float4(uv, 0, 0)).rgb);
}

float3 GetExtrapolatedSingleMieScattering(
	in AtmosphereParameters atmosphere, in float4 scattering)
{
#if SHADER_API_GLCORE
	// HACK: Somehow, on GLCore, scattering.r manages to be a NAN here whenever it's written to the texture as 0, even tho it's not read from the texture as NAN.
	if (isnan(scattering.r))
	{
		return float3(0.0, 0.0, 0.0);
	}
#endif

	if (scattering.r == 0.0) 
	{
		return float3(0.0, 0.0, 0.0);
	}
	return scattering.rgb * scattering.a / scattering.r * (atmosphere.rayleigh_scattering.r / atmosphere.mie_scattering.r) * (atmosphere.mie_scattering / atmosphere.rayleigh_scattering);
}

IrradianceSpectrum GetCombinedScattering(
	in AtmosphereParameters atmosphere,
	in ReducedScatteringTexture scattering_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground,
	out IrradianceSpectrum single_mie_scattering) 
{
	float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(atmosphere, r, mu, mu_s, nu, ray_r_mu_intersects_ground);

	float4 combined_scattering = scatteringTextureSample(scattering_texture, uvwz);

	IrradianceSpectrum scattering = IrradianceSpectrum(combined_scattering.rgb);
	single_mie_scattering = GetExtrapolatedSingleMieScattering(atmosphere, combined_scattering);

	return scattering;
}

RadianceSpectrum GetSkyRadiance(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	in ReducedScatteringTexture scattering_texture,
	Position camera, in Direction view_ray, Length shadow_length,
	in Direction sun_direction, out DimensionlessSpectrum transmittance) 
{
	// Compute the distance to the top atmosphere boundary along the view ray,
	// assuming the viewer is in space (or NaN if the view ray does not intersect
	// the atmosphere).
	Length r = length(camera);
	Length rmu = dot(camera, view_ray);
	Length distance_to_top_atmosphere_boundary = -rmu - sqrt(rmu * rmu - r * r + atmosphere.top_radius * atmosphere.top_radius);
	// If the viewer is in space and the view ray intersects the atmosphere, move
	// the viewer to the top atmosphere boundary (along the view ray):
	if (distance_to_top_atmosphere_boundary > 0.0 * m) 
	{
		camera = camera + view_ray * distance_to_top_atmosphere_boundary;
		r = atmosphere.top_radius;
		rmu += distance_to_top_atmosphere_boundary;
	}
	else if (r > atmosphere.top_radius)  // If the view ray does not intersect the atmosphere, simply return 0.
	{
		transmittance = DimensionlessSpectrum(1.0, 1.0, 1.0);
		return RadianceSpectrum(0.0, 0.0, 0.0);
	}
	// Compute the r, mu, mu_s and nu parameters needed for the texture lookups.
	Number mu = rmu / r;
	Number mu_s = dot(camera, sun_direction) / r;
	Number nu = dot(view_ray, sun_direction);
	bool ray_r_mu_intersects_ground = RayIntersectsGround(atmosphere, r, mu);

	transmittance = ray_r_mu_intersects_ground ? DimensionlessSpectrum(0.0, 0.0, 0.0) :
		GetTransmittanceToTopAtmosphereBoundary(
			atmosphere, transmittance_texture, r, mu);
	IrradianceSpectrum single_mie_scattering;
	IrradianceSpectrum scattering;
	if (shadow_length == 0.0 * m)
	{
		scattering = GetCombinedScattering(
			atmosphere, scattering_texture,
			r, mu, mu_s, nu, ray_r_mu_intersects_ground,
			single_mie_scattering);
	}
	else 
	{
		// Case of light shafts (shadow_length is the total length noted l in our
		// paper): we omit the scattering between the camera and the _point at
		// distance l, by implementing Eq. (18) of the paper (shadow_transmittance
		// is the T(x,x_s) term, scattering is the S|x_s=x+lv term).
		Length d = shadow_length;
		Length r_p =
			ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
		Number mu_p = (r * mu + d) / r_p;
		Number mu_s_p = (r * mu_s + d * nu) / r_p;

		scattering = GetCombinedScattering(
			atmosphere, scattering_texture,
			r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
			single_mie_scattering);
		DimensionlessSpectrum shadow_transmittance =
			GetTransmittance(atmosphere, transmittance_texture,
				r, mu, shadow_length, ray_r_mu_intersects_ground);
		scattering = scattering * shadow_transmittance;
		single_mie_scattering = single_mie_scattering * shadow_transmittance;
	}
	return scattering * RayleighPhaseFunction(nu) + single_mie_scattering * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
}

RadianceSpectrum GetSkyRadianceToPoint(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	in ReducedScatteringTexture scattering_texture,
	Position camera, in Direction _cameraToPoint, Length _cameraToPointDistance, Length shadow_length,
	in Direction sun_direction, out DimensionlessSpectrum transmittance)
{
	// Compute the distance to the top atmosphere boundary along the view ray,
	// assuming the viewer is in space (or NaN if the view ray does not intersect
	// the atmosphere).
	Direction view_ray = normalize(_cameraToPoint);
	Length r = length(camera);
	Length rmu = dot(camera, view_ray);
	Length distance_to_top_atmosphere_boundary = -rmu -
		sqrt(rmu * rmu - r * r + atmosphere.top_radius * atmosphere.top_radius);
	// If the viewer is in space and the view ray intersects the atmosphere, move
	// the viewer to the top atmosphere boundary (along the view ray):
	if (distance_to_top_atmosphere_boundary > 0.0 * m) {
		camera = camera + view_ray * distance_to_top_atmosphere_boundary;
		r = atmosphere.top_radius;
		rmu += distance_to_top_atmosphere_boundary;
	}

	// Compute the r, mu, mu_s and nu parameters for the first texture lookup.
	Number mu = rmu / r;
	Number mu_s = dot(camera, sun_direction) / r;
	Number nu = dot(view_ray, sun_direction);
	Length d = _cameraToPointDistance;
	bool ray_r_mu_intersects_ground = RayIntersectsGround(atmosphere, r, mu);

	transmittance = GetTransmittance(atmosphere, transmittance_texture,
		r, mu, d, ray_r_mu_intersects_ground);

	IrradianceSpectrum single_mie_scattering;
	IrradianceSpectrum scattering = GetCombinedScattering(
		atmosphere, scattering_texture,
		r, mu, mu_s, nu, ray_r_mu_intersects_ground,
		single_mie_scattering);

	// Compute the r, mu, mu_s and nu parameters for the second texture lookup.
	// If shadow_length is not 0 (case of light shafts), we want to ignore the
	// scattering along the last shadow_length meters of the view ray, which we
	// do by subtracting shadow_length from d (this way scattering_p is equal to
	// the S|x_s=x_0-lv term in Eq. (17) of our paper).
	d = max(d - shadow_length, 0.0 * m);
	Length r_p = ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
	Number mu_p = (r * mu + d) / r_p;
	Number mu_s_p = (r * mu_s + d * nu) / r_p;

	IrradianceSpectrum single_mie_scattering_p;
	IrradianceSpectrum scattering_p = GetCombinedScattering(
		atmosphere, scattering_texture,
		r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
		single_mie_scattering_p);

	// Combine the lookup results to get the scattering between camera and _point.
	DimensionlessSpectrum shadow_transmittance = transmittance;
	if (shadow_length > 0.0 * m) {
		// This is the T(x,x_s) term in Eq. (17) of our paper, for light shafts.
		shadow_transmittance = GetTransmittance(atmosphere, transmittance_texture,
			r, mu, d, ray_r_mu_intersects_ground);
	}
	scattering = scattering - shadow_transmittance * scattering_p;
	single_mie_scattering =
		single_mie_scattering - shadow_transmittance * single_mie_scattering_p;

	single_mie_scattering = GetExtrapolatedSingleMieScattering(
		atmosphere, float4(scattering, single_mie_scattering.r));

	// Hack to avoid rendering artifacts when the sun is below the horizon.
	single_mie_scattering = single_mie_scattering *
		smoothstep(Number(0.0), Number(0.01), mu_s);

	return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
		MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
}

IrradianceSpectrum GetSunAndSkyIrradiance(
	in AtmosphereParameters atmosphere,
	in TransmittanceTexture transmittance_texture,
	in IrradianceTexture irradiance_texture,
	in Position _point, in Direction normal, in Direction sun_direction,
	out IrradianceSpectrum sky_irradiance) 
{
	Length r = length(_point);
	Number mu_s = dot(_point, sun_direction) / r;

	// Indirect irradiance (approximated if the surface is not horizontal).
	sky_irradiance = GetIrradiance(atmosphere, irradiance_texture, r, mu_s) *
		(1.0 + dot(normal, _point) / r) * 0.5;

	// Direct irradiance.
	return atmosphere.solar_irradiance *
		GetTransmittanceToSun(
			atmosphere, transmittance_texture, r, mu_s) *
		smoothstep(-atmosphere.sun_angular_radius / rad,
			atmosphere.sun_angular_radius / rad,
			mu_s) *
		max(dot(normal, sun_direction), 0.0);
}
#endif