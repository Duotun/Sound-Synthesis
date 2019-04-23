using UnityEngine;
using System.Collections.Generic;

namespace PhysSound
{
    [AddComponentMenu("PhysSound/PhysSound Object")]   //add anywhere
    public class PhysSoundObject : PhysSoundBase
    {
        public PhysSoundMaterial SoundMaterial;

        public bool AutoCreateSources;
        public bool PlayClipAtPoint;
        public bool HitsTriggers;

        public AudioSource ImpactAudio;

        public List<PhysSoundAudioContainer> AudioContainers = new List<PhysSoundAudioContainer>();

        private float baseImpactVol, baseImpactPitch;
     
        private Dictionary<int, PhysSoundAudioContainer> _audioContainersDic;

        private Vector3 _prevVelocity;
        private bool _setPrevVelocity = true;

        private Vector3 _prevPosition;
        private Vector3 _kinematicVelocity;
        private Quaternion _prevRotation;
        private float _kinematicAngularVelocity;

        private int _lastFrame;

        private Rigidbody _r;
        private Rigidbody2D _r2D;

        void Start()
        {
            if (SoundMaterial == null)
                return;

            Initialize();
        }

        /// <summary>
        /// Initializes the PhysSoundObject. Use this if you adding a PhysSoundObject component to an object at runtime.
        /// </summary>
        public void Initialize()
        {
            _r = GetComponent<Rigidbody>();
            _r2D = GetComponent<Rigidbody2D>();

            if (AutoCreateSources)
            {
                baseImpactVol = ImpactAudio.volume;
                baseImpactPitch = ImpactAudio.pitch;

                _audioContainersDic = new Dictionary<int, PhysSoundAudioContainer>();
                AudioContainers = new List<PhysSoundAudioContainer>();

                foreach (PhysSoundAudioSet audSet in SoundMaterial.AudioSets)
                {
                    if (audSet.Slide == null)
                        continue;

                    PhysSoundAudioContainer audCont = new PhysSoundAudioContainer(audSet.Key);
                    audCont.SlideAudio = PhysSoundTempAudioPool.GetAudioSourceCopy(ImpactAudio, this.gameObject);

                    audCont.Initialize(SoundMaterial);
                    _audioContainersDic.Add(audCont.KeyIndex, audCont);
                    AudioContainers.Add(audCont);
                }

                ImpactAudio.loop = false;
            }
            else
            {
                if (ImpactAudio)
                {
                    ImpactAudio.loop = false;
                    baseImpactVol = ImpactAudio.volume;
                    baseImpactPitch = ImpactAudio.pitch;
                }

                if (AudioContainers.Count > 0)
                {
                    _audioContainersDic = new Dictionary<int, PhysSoundAudioContainer>();

                    foreach (PhysSoundAudioContainer audCont in AudioContainers)
                    {
                        if (!SoundMaterial.HasAudioSet(audCont.KeyIndex))
                        {
                            Debug.LogError("PhysSound Object " + gameObject.name + " has an audio container for an invalid Material Type! Select this object in the hierarchy to update its audio container list.");
                            continue;
                        }

                        audCont.Initialize(SoundMaterial);
                        _audioContainersDic.Add(audCont.KeyIndex, audCont);
                    }
                }
            }

            if (PlayClipAtPoint)
                PhysSoundTempAudioPool.Create();
            else if (ImpactAudio != null && !ImpactAudio.isActiveAndEnabled)
                ImpactAudio = PhysSoundTempAudioPool.GetAudioSourceCopy(ImpactAudio, gameObject);
        }

        void Update()
        {
            if (SoundMaterial == null)
                return;

            for (int i = 0; i < AudioContainers.Count; i++)
                AudioContainers[i].UpdateVolume();

            if (ImpactAudio && !ImpactAudio.isPlaying)
                ImpactAudio.Stop();

            _kinematicVelocity = (transform.position - _prevPosition) / Time.deltaTime;
            _prevPosition = transform.position;

            _kinematicAngularVelocity = Quaternion.Angle(_prevRotation, transform.rotation) / Time.deltaTime / 45f;
            _prevRotation = transform.rotation;        
        }

        /// <summary>
        /// Enables or Disables this script along with its associated AudioSources.
        /// </summary>
        public void SetEnabled(bool enable)
        {
            if (enable && this.enabled == false)
            {
                for (int i = 0; i < AudioContainers.Count; i++)
                {
                    AudioContainers[i].Enable();
                }

                ImpactAudio.enabled = true;
                this.enabled = true;
            }
            else if(!enable && this.enabled == true)
            {
                if (ImpactAudio)
                {
                    ImpactAudio.Stop();
                    ImpactAudio.enabled = false;
                }

                for (int i = 0; i < AudioContainers.Count; i++)
                {
                    AudioContainers[i].Disable();
                }

                this.enabled = false;
            }
        }

        /// <summary>
        /// Gets the PhysSound Material of this object.
        /// </summary>
        public override PhysSoundMaterial GetPhysSoundMaterial(Vector3 contactPoint)
        {
            return SoundMaterial;
        }

        private Vector3 TotalKinematicVelocity
        {
            get { return _kinematicVelocity + (Vector3.one * _kinematicAngularVelocity); }
        }

        #region Main Functions

        private void playImpactSound(GameObject otherObject, Vector3 relativeVelocity, Vector3 normal, Vector3 contactPoint)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || Time.frameCount == _lastFrame)
            {
                return;
            }
            
            if (ImpactAudio)
            {
                AudioClip a = SoundMaterial.GetImpactAudio(otherObject, relativeVelocity, normal, contactPoint);

                if (a)
                {
                    float pitch = baseImpactPitch * SoundMaterial.GetScaleModPitch(transform.localScale) + SoundMaterial.GetRandomPitch();
                    float vol = baseImpactVol * SoundMaterial.GetScaleModVolume(transform.localScale) * SoundMaterial.GetImpactVolume(relativeVelocity, normal);

                    if (PlayClipAtPoint)
                    {
                        PhysSoundTempAudioPool.Instance.PlayClip(a, transform.position, ImpactAudio, SoundMaterial.ScaleImpactVolume ? vol : ImpactAudio.volume, pitch);
                    }
                    else
                    {
                        ImpactAudio.pitch = pitch;
                        if (SoundMaterial.ScaleImpactVolume)
                            ImpactAudio.volume = vol;

                        ImpactAudio.clip = a;
                        ImpactAudio.Play();
                    }

                    _lastFrame = Time.frameCount;
                }
            }
        }

        private void setSlideTargetVolumes(GameObject otherObject, Vector3 relativeVelocity, Vector3 normal, Vector3 contactPoint, bool exit)
        {
            //log("Sliding! " + gameObject.name + " against " + otherObject.name + " - Relative Velocity: " + relativeVelocity + ", Normal: " + normal + ", Contact Point: " + contactPoint + ", Exit: " + exit);

            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0)
            {
                return;
            }

            PhysSoundMaterial m = null;
            PhysSoundBase b = otherObject == null ? null : otherObject.GetComponent<PhysSoundBase>();

            if (b)
            {
                //Special case for sliding against a terrain
                if (b is PhysSoundTerrain)
                {
                    PhysSoundTerrain terr = b as PhysSoundTerrain;
                    Dictionary<int, PhysSoundComposition> compDic = terr.GetComposition(contactPoint);

                    foreach (PhysSoundAudioContainer c in _audioContainersDic.Values)
                    {
                        PhysSoundComposition comp;
                        float mod = 0;

                        if (compDic.TryGetValue(c.KeyIndex, out comp))
                            mod = comp.GetAverage();

                        c.SetTargetVolumeAndPitch(this.gameObject, otherObject, relativeVelocity, normal, exit, mod);
                    }

                    return;
                }
                else
                    m = b.GetPhysSoundMaterial(contactPoint);
            }

            //General cases
            //If the other object has a PhysSound material
            if (m)
            {
                PhysSoundAudioContainer aud;

                if (_audioContainersDic.TryGetValue(m.MaterialTypeKey, out aud))
                    aud.SetTargetVolumeAndPitch(this.gameObject, otherObject, relativeVelocity, normal, exit);
                else if (!SoundMaterial.HasAudioSet(m.MaterialTypeKey) && SoundMaterial.FallbackTypeKey != -1 && _audioContainersDic.TryGetValue(SoundMaterial.FallbackTypeKey, out aud))
                    aud.SetTargetVolumeAndPitch(this.gameObject, otherObject, relativeVelocity, normal, exit);
            }
            //If it doesnt we set volumes based on the fallback setting of our material
            else
            {
                PhysSoundAudioContainer aud;

                if (SoundMaterial.FallbackTypeKey != -1 && _audioContainersDic.TryGetValue(SoundMaterial.FallbackTypeKey, out aud))
                    aud.SetTargetVolumeAndPitch(this.gameObject, otherObject, relativeVelocity, normal, exit);
            }
        }

        #endregion

        #region 3D Collision Messages

        Vector3 contactNormal, contactPoint, relativeVelocity;

        void OnCollisionEnter(Collision c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0)
                return;

            contactNormal = c.contacts[0].normal;
            contactPoint = c.contacts[0].point;
            relativeVelocity = c.relativeVelocity;

            playImpactSound(c.collider.gameObject, relativeVelocity, contactNormal, contactPoint);

            _setPrevVelocity = true;
        }


        void OnCollisionStay(Collision c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || _audioContainersDic == null)
                return;

            if (_setPrevVelocity)
            {
                _prevVelocity = _r.velocity;
                _setPrevVelocity = false;
            }

            Vector3 deltaVel = _r.velocity - _prevVelocity;

            contactNormal = c.contacts[0].normal;
            contactPoint = c.contacts[0].point;
            relativeVelocity = c.relativeVelocity;

            playImpactSound(c.collider.gameObject, deltaVel, contactNormal, contactPoint);
            setSlideTargetVolumes(c.collider.gameObject, relativeVelocity, contactNormal, contactPoint, false);

            _prevVelocity = _r.velocity;
        }

        void OnCollisionExit()
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || _audioContainersDic == null)
                return;

            setSlideTargetVolumes(null, Vector3.zero, Vector3.zero, transform.position, true);
            _setPrevVelocity = true;
        }

        #endregion

        #region 3D Trigger Messages

        void OnTriggerEnter(Collider c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || !HitsTriggers)
                return;

            playImpactSound(c.gameObject, TotalKinematicVelocity, Vector3.zero, c.transform.position);
        }

        void OnTriggerStay(Collider c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || _audioContainersDic == null || !HitsTriggers)
                return;

            setSlideTargetVolumes(c.gameObject, TotalKinematicVelocity, Vector3.zero, c.transform.position, false);
        }

        void OnTriggerExit(Collider c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || _audioContainersDic == null || !HitsTriggers)
                return;

            setSlideTargetVolumes(c.gameObject, TotalKinematicVelocity, Vector3.zero, c.transform.position, true);
        }

        #endregion

        #region 2D Collision Messages

        void OnCollisionEnter2D(Collision2D c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0)
                return;

            playImpactSound(c.collider.gameObject, c.relativeVelocity, c.contacts[0].normal, c.contacts[0].point);

            _setPrevVelocity = true;
        }

        void OnCollisionStay2D(Collision2D c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || _audioContainersDic == null)
                return;

            if (_setPrevVelocity)
            {
                _prevVelocity = _r2D.velocity;
                _setPrevVelocity = false;
            }

            Vector3 deltaVel = _r2D.velocity - (Vector2)_prevVelocity;

            playImpactSound(c.collider.gameObject, deltaVel, c.contacts[0].normal, c.contacts[0].point);
            setSlideTargetVolumes(c.collider.gameObject, c.relativeVelocity, c.contacts[0].normal, c.contacts[0].point, false);

            _prevVelocity = _r2D.velocity;
        }

        void OnCollisionExit2D(Collision2D c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || _audioContainersDic == null)
                return;

            setSlideTargetVolumes(c.collider.gameObject, c.relativeVelocity, Vector3.up, transform.position, true);

            _setPrevVelocity = true;
        }

        #endregion

        #region 2D Trigger Messages

        void OnTriggerEnter2D(Collider2D c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0)
                return;

            playImpactSound(c.gameObject, TotalKinematicVelocity, Vector3.zero, c.transform.position);
        }

        void OnTriggerStay2D(Collider2D c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || _audioContainersDic == null)
                return;

            setSlideTargetVolumes(c.gameObject, TotalKinematicVelocity, Vector3.zero, c.transform.position, false);
        }

        void OnTriggerExit2D(Collider2D c)
        {
            if (SoundMaterial == null || !this.enabled || SoundMaterial.AudioSets.Count == 0 || _audioContainersDic == null)
                return;

            setSlideTargetVolumes(c.gameObject, TotalKinematicVelocity, Vector3.zero, c.transform.position, true);
        }

        #endregion

        #region Editor

        public bool HasAudioContainer(int keyIndex)
        {
            foreach (PhysSoundAudioContainer aud in AudioContainers)
            {
                if (aud.CompareKeyIndex(keyIndex))
                    return true;
            }

            return false;
        }

        public void AddAudioContainer(int keyIndex)
        {
            AudioContainers.Add(new PhysSoundAudioContainer(keyIndex));
        }

        public void RemoveAudioContainer(int keyIndex)
        {
            for (int i = 0; i < AudioContainers.Count; i++)
            {
                if (AudioContainers[i].KeyIndex == keyIndex)
                {
                    AudioContainers.RemoveAt(i);
                    return;
                }
            }
        }

        #endregion


    }

    [System.Serializable]
    public class PhysSoundAudioContainer
    {
        public int KeyIndex;
        public AudioSource SlideAudio;

        private PhysSoundMaterial _mat;
        private float _targetVolume;
        private float _baseVol, _basePitch, _basePitchRand;

        private int _lastFrame;

        public PhysSoundAudioContainer(int k)
        {
            KeyIndex = k;
        }

        /// <summary>
        /// Initializes this Audio Container with the given AudioClip. Will do nothing if SlideAudio is not assigned.
        /// </summary>
        /// <param name="clip"></param>
        public void Initialize(PhysSoundMaterial m)
        {
            if (SlideAudio == null)
                return;

            _mat = m;

            SlideAudio.clip = _mat.GetAudioSet(KeyIndex).Slide;
            _baseVol = SlideAudio.volume;
            _basePitch = SlideAudio.pitch;
            _basePitchRand = _basePitch;
            SlideAudio.loop = true;
            SlideAudio.volume = 0;
        }

        /// <summary>
        /// Sets the target volume and pitch of the sliding sound effect based on the given object that was hit, velocity, and normal.
        /// </summary>
        public void SetTargetVolumeAndPitch(GameObject parentObject, GameObject otherObject, Vector3 relativeVelocity, Vector3 normal, bool exit, float mod = 1)
        {
            if (SlideAudio == null)
                return;

            float vol = exit || !_mat.CollideWith(otherObject) ? 0 : _mat.GetSlideVolume(relativeVelocity, normal) * _baseVol * mod;

            if (_lastFrame == Time.frameCount)
            {
                if(_targetVolume < vol)
                    _targetVolume = vol;
            }
            else
                _targetVolume = vol;

            if (!SlideAudio.isPlaying)
            {
                _basePitchRand = _basePitch * _mat.GetScaleModPitch(parentObject.transform.localScale)  + _mat.GetRandomPitch();
                SlideAudio.Play();
            }

            SlideAudio.pitch = _basePitchRand + relativeVelocity.magnitude * _mat.SlidePitchMod;

            _lastFrame = Time.frameCount;
        }

        /// <summary>
        /// Updates the associated AudioSource to match the target volume and pitch.
        /// </summary>
        public void UpdateVolume()
        {
            if (SlideAudio == null)
                return;

            SlideAudio.volume = Mathf.MoveTowards(SlideAudio.volume, _targetVolume, 0.06f);

            if (SlideAudio.volume < 0.01f)
                SlideAudio.Stop();
        }

        /// <summary>
        /// Returns true if this Audio Container's key index is the same as the given key index.
        /// </summary>
        public bool CompareKeyIndex(int k)
        {
            return k == KeyIndex;
        }

        /// <summary>
        /// Disables the associated AudioSource.
        /// </summary>
        public void Disable()
        {
            if (SlideAudio)
            {
                SlideAudio.Stop();
                SlideAudio.enabled = false;
            }
        }

        /// <summary>
        /// Enables the associated AudioSource.
        /// </summary>
        public void Enable()
        {
            if (SlideAudio)
            {
                SlideAudio.enabled = true;
            }
        }
    }
}