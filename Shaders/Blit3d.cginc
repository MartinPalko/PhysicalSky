#ifndef BLIT_3D
#define BLIT_3D

struct appdata_img3d
{
	float4 vertex : POSITION;
	half3 texcoord : TEXCOORD0;
};

struct v2f_img3d
{
	float4 pos : SV_POSITION;
	half3 uvw : TEXCOORD0;
};

v2f_img3d vert_img3d( appdata_img3d v )
{
	v2f_img3d o;
	o.pos = UnityObjectToClipPos (v.vertex);
	o.uvw = v.texcoord;
	return o;
}
#endif