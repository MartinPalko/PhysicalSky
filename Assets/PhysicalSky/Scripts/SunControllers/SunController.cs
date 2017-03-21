using UnityEngine;

namespace PhysicalSky
{
    public abstract class SunController : MonoBehaviour
    {
        protected IPhysicalSky sky;

        protected virtual void Start()
        {
            sky = GetComponent(typeof(IPhysicalSky)) as IPhysicalSky;
        }
    }
}