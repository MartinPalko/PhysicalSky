using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PhysicalSky : MonoBehaviour
{
    public AtmosphereModel atmosphereModel;

    private void Start()
    {
        atmosphereModel = GetComponent<AtmosphereModel>();
        atmosphereModel.ComputeLookupTextures();

        GetComponent<MeshRenderer>().material.SetTexture("_MainTex", atmosphereModel.ScatteringLUT);
    }

}
