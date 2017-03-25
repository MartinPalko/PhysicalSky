using UnityEngine;

namespace PhysicalSky
{
    [ExecuteInEditMode]
    public abstract class SunController : MonoBehaviour
    {
        protected IPhysicalSky sky;

        protected virtual void Start()
        {
            sky = GetComponent(typeof(IPhysicalSky)) as IPhysicalSky;
        }
    }
}