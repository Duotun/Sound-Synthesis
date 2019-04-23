using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;
using PhysSound;
public class Setoff : MonoBehaviour {
    public GameObject[] Objects;

    public Transform DropLocation;
    public float RandomForce;

    public Material soft;
    public Material hard;
    GameObject ground;
    public PhysSoundMaterial SoundMat_hard;
    public PhysSoundMaterial SoundMat_soft;
    void Start()
    {
        foreach (GameObject g in Objects)
            g.GetComponent<Rigidbody>().maxAngularVelocity = 1000;
        ground = GameObject.Find("Floor");
        drop(-1);
    }

    void Update()
    {
        if (Input.GetKeyDown(KeyCode.D))
        {
            drop(-1);
        }

        if (Input.GetKeyDown(KeyCode.Alpha1))
            drop(0);
        if (Input.GetKeyDown(KeyCode.Alpha2))
            drop(1);
        if (Input.GetKeyDown(KeyCode.Alpha3))
            drop(2);
        if (Input.GetKeyDown(KeyCode.Alpha4))
            drop(3);
        if (Input.GetKeyDown(KeyCode.Alpha5))
            drop(4);
        if(Input.GetKeyDown(KeyCode.S))
        {
            change_to_soft();
        }
        if(Input.GetKeyDown(KeyCode.F))
        {
            change_to_hard();
        }
    }

    void change_to_soft()  //change materials
    {
        ground.GetComponent<PhysSoundObject>().SoundMaterial = SoundMat_soft;
        ground.GetComponent<MeshRenderer>().material = soft;
    }
    void change_to_hard()
    {
        ground.GetComponent<PhysSoundObject>().SoundMaterial = SoundMat_hard;
        ground.GetComponent<MeshRenderer>().material = hard;

    }
    void drop(int obj)
    {
        for (int i = 0; i < Objects.Length; i++)
        {
            if (obj == -1 || i == obj)
            {
                GameObject g = Objects[i];
                g.transform.position = DropLocation.position + Random.insideUnitSphere * 1.5f;
                g.GetComponent<Rigidbody>().velocity = new Vector3(Random.Range(-RandomForce, RandomForce), Random.Range(-RandomForce, RandomForce), 0);
                g.GetComponent<Rigidbody>().angularVelocity = Random.insideUnitSphere * RandomForce;
            }
        }
    }

    void OnGUI()
    {
        GUILayout.BeginArea(new Rect(10, 10, Screen.width - 20, 50));
        GUILayout.BeginHorizontal();
        GUILayout.FlexibleSpace();
        GUILayout.Box("Sound Testing Demo");   //SceneManager.GetActiveScene().name
        GUILayout.Box("Press 'D' to drop objects.");
        GUILayout.Box("Press '1', '2', '3' to drop specific objects.");
        GUILayout.Box("Press 'S', 'F' to change ground materials.");
        //GUILayout.Box("Current Object: " + Target.name);
        //GUILayout.Box("'1' '2' '3' or '4' to load different scenes.");
        GUILayout.FlexibleSpace();
        GUILayout.EndHorizontal();
        GUILayout.EndArea();
    }
}
