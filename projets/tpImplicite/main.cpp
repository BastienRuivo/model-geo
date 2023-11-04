
//! \file tuto7_camera.cpp reprise de tuto7.cpp mais en derivant AppCamera, avec gestion automatique d'une camera.


#include "wavefront.h"
#include "uniforms.h"
#include "texture.h"

#include "orbiter.h"
#include "draw.h"        
#include "app_camera.h"        // classe Application a deriver

#include "implicit.h"
#include <omp.h>
#include <chrono>

#include "Bezier.h"
#include "Deformation.h"




BOOST_CLASS_EXPORT(Point)
BOOST_CLASS_EXPORT(Vector)
BOOST_CLASS_EXPORT(vec2)
BOOST_CLASS_EXPORT(vec3)
BOOST_CLASS_EXPORT(vec4)
BOOST_CLASS_EXPORT(Box)


BOOST_CLASS_EXPORT(Node)
BOOST_CLASS_EXPORT(BinaryOperator)
BOOST_CLASS_EXPORT(Tree)
BOOST_CLASS_EXPORT(Union)
BOOST_CLASS_EXPORT(UnionSmooth)
BOOST_CLASS_EXPORT(Difference)
BOOST_CLASS_EXPORT(DifferenceSmooth)
BOOST_CLASS_EXPORT(Intersection)
BOOST_CLASS_EXPORT(IntersectionSmooth)
BOOST_CLASS_EXPORT(Sphere)
BOOST_CLASS_EXPORT(Cube)
BOOST_CLASS_EXPORT(Tor)
BOOST_CLASS_EXPORT(Cylinder)
BOOST_CLASS_EXPORT(Capsule)

enum ray_algo {
    RAY_ALGO_SPHERE_TRACING,
    RAY_ALGO_SPHERE_TRACING_BBOX,
    RAY_ALGO_RAY_MARCHING,
    RAY_ALGO_RAY_MARCHING_BBOX
};

enum primitive {
    PRIMITIVE_BOX,
    PRIMITIVE_SPHERE,
    PRIMITIVE_CAPSULE,
    PRIMITIVE_TOR
};

enum binary_operator {
    BINARY_OPERATOR_UNION,
    BINARY_OPERATOR_UNION_SMOOTH,
    BINARY_OPERATOR_INTERSECTION,
    BINARY_OPERATOR_INTERSECTION_SMOOTH,
    BINARY_OPERATOR_DIFFERENCE,
    BINARY_OPERATOR_DIFFERENCE_SMOOTH
};




const float speed= 0.075f;        // facteur d'avancement automatique

// utilitaire. creation d'une grille / repere.
Mesh make_grid( const int n= 10 )
{
    Mesh grid= Mesh(GL_LINES);
    
    // grille
    grid.color(White());
    for(int x= 0; x < n; x++)
    {
        float px= float(x) - float(n)/2 + .5f;
        grid.vertex(Point(px, 0, - float(n)/2 + .5f)); 
        grid.vertex(Point(px, 0, float(n)/2 - .5f));
    }

    for(int z= 0; z < n; z++)
    {
        float pz= float(z) - float(n)/2 + .5f;
        grid.vertex(Point(- float(n)/2 + .5f, 0, pz)); 
        grid.vertex(Point(float(n)/2 - .5f, 0, pz)); 
    }
    
    // axes XYZ
    // grid.color(Red());
    // grid.vertex(Point(0, .1, 0));
    // grid.vertex(Point(1, .1, 0));
    
    // grid.color(Green());
    // grid.vertex(Point(0, .1, 0));
    // grid.vertex(Point(0, 1, 0));
    
    // grid.color(Blue());
    // grid.vertex(Point(0, .1, 0));
    // grid.vertex(Point(0, .1, 1));
    
    glLineWidth(2);
    
    return grid;
}


class TP : public AppCamera
{
    void makeTree(Node * root) {
        if(groot != nullptr) delete groot;
        groot = new Tree(root);
        updateTree();
    }

    void makeTree(const std::string & filename) {
        if(groot != nullptr) delete groot;
        groot = new Tree(filename);
        updateTree();
    }
    void updateTree() {
        auto start = std::chrono::high_resolution_clock::now();
        nbPrimitive = groot->countPrimitives();
        m_implicit= Mesh(GL_TRIANGLES);
        groot->Polygonize(resolution, m_implicit, Space);

        auto end = std::chrono::high_resolution_clock::now();
        lastComputeTime = std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count();
        lastVerticesCount = m_implicit.vertex_count();
        lastTrianglesCount = m_implicit.triangle_count();
    }
    void makePortalCube() {

        Node * t = new Cube(Vector(0.5, 0.5, 0.5), 1.0);
        t = new Difference(t, new Sphere(Vector(0.5, 0.5, 0.5), 0.725));
        t = new Union(t, new Cube(Vector(0.5, 0.5, 0.5), 0.97));
        t = new Difference(t, new Sphere(Vector(0.5, 0.5, 0.5), 0.6));
        t = new Union(t, new Cube(Vector(0.5, 0.5, 0.5), 0.94));
        t = new Union(t, new Sphere(Vector(0.5, 0.5, 0.5), 0.5));

        makeTree(t);
    }
public:
    // constructeur : donner les dimensions de l'image, et eventuellement la version d'openGL.
    TP( ) : AppCamera(1900, 1000) {}
    
    // creation des objets de l'application
    int init( )
    {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();

        io = ImGui::GetIO();(void)io;
        ImGui_ImplSdlGL3_Init(m_window);

        ImGui::StyleColorsDark();

        // decrire un repere / grille 
        m_repere= make_grid(10);
        
        // charge un objet
        m_cube= read_mesh("data/cube.obj");
        
        shaderProg = read_program("projets/tp/shader.glsl");
        program_print_errors(shaderProg);

        //makeTree(std::string("data/") + std::string(filename) + std::string(".txt"));
        makeTree(new Cube(Vector(0, 0, 0), 1));


        

        std::vector<std::vector<Point>> pointVase {
            {Point(0, 0, 0), Point(0.1, 0.0, 0.0), Point(0.15, 0.0, 0.0), Point(0.2, 0.1, 0.0), Point(0.3, 0.6, 0.0), Point(0.2, 0.7, 0.0), Point(0.1, 0.7, 0.0), Point(0.1, 0.9, 0.0), Point(0.25, 0.95, 0.0)}
        };
        vase = Bezier(pointVase);
        Revolution vaseRev(vase, Vector(0, 1, 0));

        FastNoiseLite noise;
        noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
        noise.SetFrequency(2.f);
        noise.SetSeed(1337);
        std::vector<std::vector<Point>> pointSurface;
        int resnoise = 32;
        for (size_t i = 0; i < resnoise; i++)
        {
            pointSurface.push_back(std::vector<Point>());
            for (size_t j = 0; j < resnoise; j++)
            {
                float x = i / (float)resnoise;
                float z = j / (float)resnoise;
                float y = (noise.GetNoise(x, z) + 1) * 0.5;
                Point p = Point(x, y, z);
                pointSurface[i].push_back(Point(x, y, z));
            }
            
        }

        bezier = Bezier(pointSurface);
        int res = 4;
        // for (size_t i = 0; i < 8; i++)
        // {
            
        //     std::cout<<"Polygonizing res = " << res * (i + 1)<<std::endl;
        //     m_Beziers.push_back(bezier.polygonize((i+1) * res));
        //     m_Revolutions.push_back(vaseRev.polygonize((i+1) * res));
        // }
        // std::cout<<"Polygonizing done"<<std::endl;

        for (size_t i = 0; i < 4; i++)
        {
            m_Beziers.push_back(GlobalDeformation::TorsionHelicoidal::Warp(m_implicit, 1, 1.0 * (i + 1)));
        }
        
        
        

        srand(time(NULL));


        std::cout<< "Getting lines mesh"<<std::endl;
        lines_bezier.clear();
        lines_revolution.clear();
        bezier.getLinesMesh(lines_bezier);
        vase.getLinesMesh(lines_revolution);
        std::cout<< "Getting lines mesh done"<<std::endl;


        
        // etat openGL par defaut
        glClearColor(0.2f, 0.2f, 0.2f, 1.f);        // couleur par defaut de la fenetre
        
        glClearDepth(1.f);                          // profondeur par defaut
        glDepthFunc(GL_LESS);                       // ztest, conserver l'intersection la plus proche de la camera
        glEnable(GL_DEPTH_TEST);                    // activer le ztest

        camera().lookat(Point(-5, -5, -5), Point(5, 5, 5)); // se positionner et regarder le point 0,0,0

        std::cout<<"INIT DONE ---- STARTING "<<std::endl;

        return 0;   // pas d'erreur, sinon renvoyer -1
    }
    
    // destruction des objets de l'application
    int quit( )
    {
        m_repere.release();
        return 0;   // pas d'erreur
    }
    void doImplicitUi() {
        ImGui::Begin("UI");
        ImGui::SeparatorText("Formes implicites");
        if(ImGui::Button("Box")) {
            Node * b = new Cube();
            makeTree(b);
        }
        ImGui::SameLine();
        if(ImGui::Button("Sphere")) {
            Node * sph = new Sphere(Vector(0, 0, 0), 0.5);
            makeTree(sph);
        }
        ImGui::SameLine();
        if(ImGui::Button("Capsule")) {
            Node * cap = new Capsule(Vector(0, 0, 0), 0.25, 0.5);
            makeTree(cap);
        }
        ImGui::SameLine();
        if(ImGui::Button("Tor")) {
            Node * tor = new Tor(Vector(0, 0, 0), 0.25, 0.5);
            makeTree(tor);
        }
        ImGui::SameLine();
        if(ImGui::Button("Portal cube")) {
            makePortalCube();
        }
        ImGui::SeparatorText("Marching Cube Parameters");
        ImGui::InputInt("Resolution", &resolution);
        ImGui::SameLine();
        if(ImGui::Button("remesh")) {
            updateTree();
        }
        ImGui::InputFloat3("Min", &(Space.a).x);
        ImGui::InputFloat3("Max", &(Space.b).x);
        ImGui::SeparatorText("Ray Tracing");
        ImGui::Combo("Algorithm", &algo, "Sphere Tracing\0Sphere Tracing BBOX\0Ray Marching\0Ray Marching BBOX\0", 4);
        
        ImGui::SeparatorText("Specs");
        ImGui::InputInt("Resolution", &resolution);
        ImGui::InputFloat3("Min", &(Space.a).x);
        ImGui::InputFloat3("Max", &(Space.b).x);
        ImGui::InputFloat("Time Rendering (s)", &lastComputeTime, 0, 0, "%.3f", ImGuiInputTextFlags_ReadOnly);
        ImGui::InputFloat("Time Intersecting (ms, 10 ray)", &lastComputeRayTime, 0, 0, "%.3f", ImGuiInputTextFlags_ReadOnly);
        ImGui::InputInt("Triangles", &lastTrianglesCount, 0, 0, ImGuiInputTextFlags_ReadOnly);
        ImGui::InputInt("Vertices", &lastVerticesCount, 0, 0, ImGuiInputTextFlags_ReadOnly);
        ImGui::InputInt("Primitives", &nbPrimitive, 0, 0, ImGuiInputTextFlags_ReadOnly);
        ImGui::SeparatorText("Scene");
        if(ImGui::Checkbox("Wireframe", &wireframe)) {
            if(wireframe) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            } else {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            }
        }
        ImGui::End();
    }   
    void doRaycastUi() {
        ImGui::Begin("Tools");
        ImGui::SeparatorText("Binary Operator");
        ImGui::Combo("Binary Operator", &binaryOperator, "Union\0Union Smooth\0Intersection\0Intersection Smooth\0Difference\0Difference Smooth\0", 6);
        if(binaryOperator == BINARY_OPERATOR_UNION_SMOOTH || binaryOperator == BINARY_OPERATOR_INTERSECTION_SMOOTH || binaryOperator == BINARY_OPERATOR_DIFFERENCE_SMOOTH) {
            ImGui::SliderFloat("Smooth Radius", &smoothRadius, 0.01f, 1.f);
        } 
        
        ImGui::SeparatorText("Primitive");
        ImGui::Combo("Primitive", &primitive, "Box\0Sphere\0Capsule\0Tor\0", 3);

        if(primitive == PRIMITIVE_BOX) {
            ImGui::InputFloat("Size", &o1);
        } else if(primitive == PRIMITIVE_SPHERE) {
            ImGui::InputFloat("Radius", &o1);
        } else if(primitive == PRIMITIVE_CAPSULE) {
            ImGui::InputFloat("Height", &o2);
            ImGui::InputFloat("Radius", &o1);
        } else if(primitive == PRIMITIVE_TOR) {
            ImGui::InputFloat("Radius Center", &o1);
            ImGui::InputFloat("Radius Tube", &o2);
        }
        ImGui::Checkbox("Auto clear ray", &autoClearRay);
        ImGui::Checkbox("Keep only addition", &keepOnlyAddtion);
        ImGui::Checkbox("Show Bbox", &showBbox);
        if(ImGui::Button("Reset Rays")) {
            m_rays.clear();
        }

        ImGui::SameLine();
        ImGui::SliderInt("Nb Rays", &nbRay, 1, 200);
        ImGui::SliderFloat("Dispersion", &dispersion, 0.01f, 1.f);

        ImGui::InputText("Filename", filename, 256);

        if(ImGui::Button("Save SDF")) {
            groot->save(std::string("data/") + std::string(filename) + std::string(".txt"));
        }

        ImGui::SameLine();

        if(ImGui::Button("Load SDF")) {
            std::cout<<"test"<<std::endl;
            groot->load(std::string("data/") + std::string(filename) + std::string(".txt"));
            updateTree();
        }

        ImGui::SameLine();

        if(ImGui::Button("Save OBJ")) {
            write_mesh(m_implicit, (std::string("./data/") + std::string(filename) + std::string(".obj")).c_str());
        }

        ImGui::End();
    }  
    void doBezierUi() {
        ImGui::Begin("Bezier");
        ImGui::Checkbox("Show Lines", &drawLines);
        ImGui::SeparatorText("Bezier");
        //ImGui::Text("Number of curves : %d", bezier.getNbCurves());
        //ImGui::Text("Number of points per curves : %d", bezier.getNbPoints());
        for (size_t i = 0; i < m_Beziers.size(); i++)
        {
            ImGui::Text("Fréquence : %d",  (i + 1));
            // ImGui::Text("Vertices : %d", m_Revolutions[i].vertex_count());
            // ImGui::Text("Triangles : %d", m_Revolutions[i].triangle_count());
        }
        

        ImGui::End();
    }
    void doUI() {
        ImGui_ImplSdlGL3_NewFrame(m_window);

        //doImplicitUi();
        //doRaycastUi();
        doBezierUi();
        ImGui::Render();
        ImGui_ImplSdlGL3_RenderDrawData(ImGui::GetDrawData());
    }

    float getRand(float range) {
        // return a float between -range and range
        return ((float)rand() / (float)(RAND_MAX)) * range * 2 - range;
    }
    
    void doKeyboard() {
        if (io.WantCaptureKeyboard) return;
        int mx, my;
        unsigned int mb = SDL_GetRelativeMouseState(&mx, &my);
        int mox, moy;
        SDL_GetMouseState(&mox, &moy);

        bool ctrlKey = key_state(SDLK_LCTRL) || key_state(SDLK_RCTRL);

        if(ctrlKey && mb && SDL_BUTTON(1)) {
            if(launched) return;
            Vector ray_nds = Vector((2.f * mox) / window_width() - 1.f, 1.f - (2.f * moy) / window_height(), 1.f);
            vec4 ray_eye = (camera().projection().inverse())(vec4(ray_nds.x, ray_nds.y, -1.0, 1.0));
            ray_eye = vec4(ray_eye.x, ray_eye.y, -1.0, 0.0);
            mPos = camera().position() + normalize((camera().view().inverse())(ray_eye)) * 5.0;

            Vector mDir = normalize((camera().view().inverse())(ray_eye));

            // Mesh ray = Mesh(GL_LINES);
            // ray.vertex(mPos - mDir * 30);

            launched = true;
            bool touched = false;
            Vector r;

            std::vector<Node*> primitives;

            if(autoClearRay) {
                m_rays.clear();
            }

            // Mesh ray = Mesh(GL_LINES);
            // ray.vertex(mPos);
            // ray.vertex(mPos + mDir * 30);
            // ray.color(Blue());
            // ray.color(Blue());
            // m_rays.push_back(ray);
            auto start = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < nbRay; i++)
            {
                // create random direction around mDir in a cone of 30 degrees
                Vector dir = mDir;
                if(nbRay != 1) {
                    dir.x = dir.x + getRand(dispersion);
                    dir.y = dir.y + getRand(dispersion);
                    dir.z = dir.z + getRand(dispersion);
                }

                Mesh ray = Mesh(GL_LINES);
                ray.vertex(mPos);

                dir = normalize(dir);
                if((algo ==  ray_algo::RAY_ALGO_SPHERE_TRACING && groot->intersectSphereTracing(mPos - dir * 30, dir, r)) ||
                (algo ==  ray_algo::RAY_ALGO_SPHERE_TRACING_BBOX && groot->intersectSphereTracingBBOX(mPos - dir * 30, dir, r)) ||
                (algo ==  ray_algo::RAY_ALGO_RAY_MARCHING && groot->intersectRayMarching(mPos - dir * 30, dir, r)) ||
                (algo ==  ray_algo::RAY_ALGO_RAY_MARCHING_BBOX && groot->intersectRayMarchingBBOX(mPos - dir * 30, dir, r)
                )) {
                    ray.vertex(r);
                    ray.color(Green());
                    ray.color(Green());
                    touched = true;
                    Node * n;
                    if(primitive == PRIMITIVE_BOX) {
                        n = new Cube(r, o1);
                    } else if(primitive == PRIMITIVE_SPHERE) {
                        n = new Sphere(r, o1);
                    } else if(primitive == PRIMITIVE_CAPSULE) {
                        n = new Capsule(r, o1, o2);
                    } else if(primitive == PRIMITIVE_TOR) {
                        n = new Tor(r, o2, o1);
                    }
                    primitives.push_back(n);
                } else {
                    ray.vertex(mPos + dir * 30);
                    ray.color(Red());
                    ray.color(Red());
                }

                m_rays.push_back(ray);
            }

            auto end = std::chrono::high_resolution_clock::now();
            lastComputeRayTime = (std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() * 1000);

            if(touched) {
                for (size_t i = 0; i < primitives.size(); i++)
                {
                     if(binaryOperator == BINARY_OPERATOR_UNION) 
                        groot->root = new Union(groot->root, primitives[i]);
                    else if(binaryOperator == BINARY_OPERATOR_UNION_SMOOTH)
                        groot->root = new UnionSmooth(groot->root, primitives[i], smoothRadius);
                    else if(binaryOperator == BINARY_OPERATOR_INTERSECTION)
                        groot->root = new Intersection(groot->root, primitives[i]);
                    else if(binaryOperator == BINARY_OPERATOR_INTERSECTION_SMOOTH)
                        groot->root = new IntersectionSmooth(groot->root, primitives[i], smoothRadius);
                    else if(binaryOperator == BINARY_OPERATOR_DIFFERENCE)
                        groot->root = new Difference(groot->root, primitives[i]);
                    else if(binaryOperator == BINARY_OPERATOR_DIFFERENCE_SMOOTH)
                        groot->root = new DifferenceSmooth(groot->root, primitives[i], smoothRadius);
                }
                updateTree();
            }

        } else {
            launched = false;
        }
    }

    void setUniforms(const Transform & model) {
        Transform view = camera().view();
        Transform projection = camera().projection();

    }

    // dessiner une nouvelle image
    int render( )
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        //set to wireframe
        // Dessinne le repere et la lumiere
        //draw(m_repere, /* model */ Identity(), camera());

        //draw(m_implicit, /* model */ Identity(), camera());
        //draw(m_implicit_deformed, Translation(Vector(m_Beziers.size() + 1.5, 0, -1)), camera());


        
        for (int i = 0; i < m_Beziers.size(); i++)
        {
            float x = m_Beziers.size() * 2 - i * 2;
            draw(m_Beziers[i], /* model */ Translation(Vector(x, 0, 0)), camera());
            if(drawLines && i == m_Beziers.size() - 1) {
                for (auto & line : lines_bezier) {
                    draw(line, /* model */ Translation(Vector(x, 0, 0)), camera());
                }
            }
        }

        for (int i = 0; i < m_Revolutions.size(); i++)
        {
            float x = m_Revolutions.size() * 2 - i * 2;
            draw(m_Revolutions[i], /* model */ Translation(Vector(x, 0, 2.5)), camera());
            if(drawLines && i == m_Revolutions.size() - 1) {
                for (auto & line : lines_revolution) {
                    draw(line, /* model */ Translation(Vector(x, 0, 2.5)), camera());
                }
            }
        }
        

        for (auto & ray : m_rays) {
            draw(ray, /* model */ Identity(), camera());
        }

        // if(showBbox && groot != nullptr && groot->root != nullptr) {
        //     std::vector<Mesh> meshes;
        //     meshes.push_back(groot->SDFBox.ToMesh());

        //     for (auto & mesh : meshes) {
        //         draw(mesh, /* model */ Identity(), camera());
        //     }
        // }


        doUI();
        doKeyboard();

        return 1;
    }

protected:
    Mesh m_cube;
    Mesh m_repere;
    Mesh m_implicit;
    Mesh m_implicit_deformed;
    std::vector<Mesh> m_rays;
    std::vector<Mesh> lines_bezier;
    std::vector<Mesh> lines_revolution;
    std::vector<Mesh> m_Beziers;
    std::vector<Mesh> m_Revolutions;

    float lastComputeRayTime = 0;

    float lastComputeTime = 0;
    int lastTrianglesCount = 0;
    int lastVerticesCount = 0;

    float smoothRadius = 1.f;

    bool wireframe = false;

    int primitive = PRIMITIVE_SPHERE;
    int binaryOperator = BINARY_OPERATOR_DIFFERENCE;

    int algo = RAY_ALGO_SPHERE_TRACING;

    float o1 = 0.1f, o2 = 0.1f, o3 = 0.1f;
    Vector aa, bb;

    bool showBbox = false;

    char filename[256] = "Bonzai";

    bool launched = false;

    bool autoClearRay = true;

    int nbRay = 1;

    float dispersion = 0.1f;

    int nbPrimitive = 1;

    Bezier bezier;

    Bezier vase;

    bool drawLines = true;


    

    Point mPos;

    bool keepOnlyAddtion = false;

    Tree * groot;

    int resolution = 72;
    Box Space = Box(Vector(-1.5, -1, -1.5), Vector(1.5, 2.5, 1.5));
    GLuint shaderProg;
    ImGuiIO io;
};


int main( int argc, char **argv )
{
    // il ne reste plus qu'a creer un objet application et la lancer 

    


    TP tp;
    tp.run();

    ImGui_ImplSdlGL3_Shutdown();
    ImGui::DestroyContext();
    
    return 0;
}
