
//! \file tuto7_camera.cpp reprise de tuto7.cpp mais en derivant AppCamera, avec gestion automatique d'une camera.


#include "wavefront.h"
#include "uniforms.h"
#include "texture.h"

#include "orbiter.h"
#include "draw.h"        
#include "app_camera.h"        // classe Application a deriver




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
    grid.color(Red());
    grid.vertex(Point(0, .1, 0));
    grid.vertex(Point(1, .1, 0));
    
    grid.color(Green());
    grid.vertex(Point(0, .1, 0));
    grid.vertex(Point(0, 1, 0));
    
    grid.color(Blue());
    grid.vertex(Point(0, .1, 0));
    grid.vertex(Point(0, .1, 1));
    
    glLineWidth(2);
    
    return grid;
}

class TP : public AppCamera
{
public:
    struct LocalMesh
    {
        GLuint VAO;
        GLuint vertex_buffer;
        GLuint uv_buffer;
        GLuint normal_buffer;
        GLuint color_buffer;
        GLuint material_buffer;
        int vertexCount;

        LocalMesh() : VAO(0), vertex_buffer(0), normal_buffer(0), material_buffer(0), vertexCount(0) {}

        void init(const Mesh & m) {
            if(!(m.vertex_buffer_size())) {
                printf("LocalMesh::init error: mesh has no vertex buffer\n");
                return;
            }

            
            glGenBuffers(1, &vertex_buffer);
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
            glBufferData(GL_ARRAY_BUFFER, m.vertex_buffer_size(), m.vertex_buffer(), GL_STATIC_DRAW); // Utilisation par draw, sans modifications

            

            if(m.has_texcoord()) {
                glGenBuffers(1, &uv_buffer);
                glBindBuffer(GL_ARRAY_BUFFER, uv_buffer);
                glBufferData(GL_ARRAY_BUFFER, m.texcoord_buffer_size(), m.texcoord_buffer(), GL_STATIC_DRAW);

                
            }

             // On recommence avec les normales

            if(m.has_normal()) {
                glGenBuffers(1, &normal_buffer);
                glBindBuffer(GL_ARRAY_BUFFER, normal_buffer);
                glBufferData(GL_ARRAY_BUFFER, m.normal_buffer_size(), m.normal_buffer(), GL_STATIC_DRAW);

                
            }
            
            // On recommence avec le material index

            if(m.has_material_index()) {
                std::vector<uint> buffer(m.vertex_count());
                for(int triangle_id= 0; triangle_id < m.triangle_count(); triangle_id++)
                {
                    int material_id = m.triangle_material_index(triangle_id);
                    unsigned a= triangle_id*3;
                    unsigned b= triangle_id*3 +1;
                    unsigned c= triangle_id*3 +2;
                    
                    buffer[a]= material_id;
                    buffer[b]= material_id;
                    buffer[c]= material_id;
                }
                
                glGenBuffers(1, &material_buffer);
                glBindBuffer(GL_ARRAY_BUFFER, material_buffer);
                glBufferData(GL_ARRAY_BUFFER, buffer.size() * sizeof(uint), buffer.data(), GL_STATIC_DRAW);
            }

            glGenVertexArrays(1, &VAO);
            glBindVertexArray(VAO);

            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
            glEnableVertexAttribArray(0);

            glBindBuffer(GL_ARRAY_BUFFER, uv_buffer);
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
            glEnableVertexAttribArray(1);

            glBindBuffer(GL_ARRAY_BUFFER, normal_buffer);
            glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);
            glEnableVertexAttribArray(2);

            glBindBuffer(GL_ARRAY_BUFFER, material_buffer);
            glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, 0, 0);
            glEnableVertexAttribArray(4);



            
            vertexCount = m.vertex_count();

            // etat openGL par defaut
            glClearColor(0.2f, 0.2f, 0.2f, 1.f); // couleur par defaut de la fenetre

            glClearDepth(1.f);       // profondeur par defaut
            glDepthFunc(GL_LESS);    // ztest, conserver l'intersection la plus proche de la camera
            glEnable(GL_DEPTH_TEST); // activer le ztest

            
        }

        void destroy() {
            glDeleteVertexArrays(1, &VAO);
            glDeleteBuffers(1, &vertex_buffer);
            glDeleteBuffers(1, &uv_buffer);
            glDeleteBuffers(1, &normal_buffer);
            glDeleteBuffers(1, &material_buffer);
        }
    };
    
    // constructeur : donner les dimensions de l'image, et eventuellement la version d'openGL.
    TP( ) : AppCamera(1900, 1000) {}
    
    // creation des objets de l'application
    int init( )
    {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();

        io = ImGui::GetIO();(void)io;
        ImGui_ImplSdlGL3_Init(m_window);
        //glewInit();

        ImGui::StyleColorsDark();

        // decrire un repere / grille 
        m_repere= make_grid(10);
        
        // charge un objet
        m_cube= read_mesh("data/cube.obj");
        
        m_robot = read_mesh("data/robot.obj");

        

        base_texture = read_texture(0, "data/terrain/Clipboard01_texture.png");

        if(m_robot.vertex_count() == 0)
        {
            printf("le fichier robot.obj n'existe pas, ou est vide\n");
            return -1;
        }

        if(m_robot.materials().count() == 0)
        {
            printf("le fichier robot.obj ne contient pas de materiaux\n");
            return -1;
        }
        
        m_materials.resize(8, White());

        m_object.init(m_robot);

        for(int i= 0; i < m_robot.materials().count(); i++) {
            m_materials[i] = m_robot.materials().material(i).diffuse;
            std::cout<<m_materials[i].r << " " << m_materials[i].g << " " << m_materials[i].b << std::endl;
        }

        std::cout<<"Robot has "<< m_robot.materials().count() <<std::endl;
        
        shaderProg = read_program("projets/tp/shader.glsl");
        program_print_errors(shaderProg);
        
        // etat openGL par defaut
        glClearColor(0.2f, 0.2f, 0.2f, 1.f);        // couleur par defaut de la fenetre
        
        glClearDepth(1.f);                          // profondeur par defaut
        glDepthFunc(GL_LESS);                       // ztest, conserver l'intersection la plus proche de la camera
        glEnable(GL_DEPTH_TEST);                    // activer le ztest

        return 0;   // pas d'erreur, sinon renvoyer -1
    }
    
    // destruction des objets de l'application
    int quit( )
    {
        m_robot.release();
        m_repere.release();
        m_object.destroy();
        return 0;   // pas d'erreur
    }

    void doUI() {
        ImGui_ImplSdlGL3_NewFrame(m_window);
        //ImGui::NewFrame();

        // render your GUI
        SDL_Event event;
        SDL_PollEvent(&event);
        
        ImGui_ImplSdlGL3_ProcessEvent(&event);
        std::cout<<io.KeyCtrl<<std::endl;
        ImGui::Begin("UI");


        ImGui::SliderFloat3("Translation", &t.x, -10, 10);
        ImGui::SliderFloat3("Scale", &s.x, -2, 2);
        ImGui::SliderAngle("Rotation X", &r.x, -180, 180);
        ImGui::SliderAngle("Rotation Y", &r.y, -180, 180);
        ImGui::SliderAngle("Rotation Z", &r.z, -180, 180);
        ImGui::ColorEdit3("Light Color", &lightColor.r);
        ImGui::SliderFloat("Light Intensity", &lightIntensity, 0, 20);
        ImGui::SliderFloat("Shininess", &shininess, 0, 20);
        

        if(ImGui::Button("Reset")) {
            t = vec3(0, 0, 0);
            s = vec3(0.1, 0.1, 0.1);
            r = vec3(0, 0, 0);
            lightColor = Color(1, 1, 1);
            lightIntensity = 10.0f;
        }

        ImGui::SameLine();

        if(ImGui::Checkbox("Wireframe", &isWireframe)) {
            if(isWireframe)
                glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
            else
                glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        }


        ImGui::End();

        ImGui::Render();
        ImGui_ImplSdlGL3_RenderDrawData(ImGui::GetDrawData());
    }
    
    void doKeyboard() {
        if (io.WantCaptureKeyboard) return;
        if(key_state(' ')) light.y += speed;
        if(key_state(SDLK_LSHIFT)) light.y -= speed;
        if(key_state('z')) light.z -= speed;
        if(key_state('s')) light.z += speed;
        if(key_state('d')) light.x += speed;
        if(key_state('q')) light.x -= speed;
        if(key_state('a')) light = vec3(0, 0, 2.5);
    }

    void setUniforms(const Transform & model) {
        Transform view = camera().view();
        Transform projection = camera().projection(window_width(), window_height(), 45);
        
        program_uniform(shaderProg, "mvpMatrix", projection * view * model);
        program_uniform(shaderProg, "modelMatrix", model);
        program_uniform(shaderProg, "source", light);
        program_uniform(shaderProg, "camera", camera().position());
        program_uniform(shaderProg, "lightColor", lightColor);
        program_uniform(shaderProg, "materials", m_materials);
        program_uniform(shaderProg, "shininess", shininess);
        program_uniform(shaderProg, "lightPower", lightIntensity);

        program_use_texture(shaderProg, "base_texture", 0, base_texture);
    }

    // dessiner une nouvelle image
    int render( )
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Dessinne le repere et la lumiere
        draw(m_repere, /* model */ Identity(), camera());
        draw(m_cube, Translation(light) * Scale(0.1, 0.1, 0.1), camera());

        // Gere les inputs clavier et souris

        // Dessine le robot
        Transform Model = 
            Translation(t.x, t.y, t.z) * 
            RotationX(r.x * 180 / M_PI) * 
            RotationY(r.y * 180 / M_PI) * 
            RotationZ(r.z * 180 / M_PI) * 
            Scale(s.x, s.y, s.z);
        
        glUseProgram(shaderProg);
        setUniforms(Model);
        glBindVertexArray(m_object.VAO);
        glDrawArrays(GL_TRIANGLES, 0, m_object.vertexCount);

        //m_robot.draw(shaderProg, true, false, true, true, true);


        doUI();
        doKeyboard();

        return 1;
    }

protected:
    Mesh m_cube;
    Mesh m_repere;

    LocalMesh m_object;

    Mesh m_robot;
    
    GLuint shaderProg;
    GLuint base_texture;

    bool isWireframe = false;

    float k = 0.8f;
    float shininess = 5.0f;

    float lightIntensity = 13.0f;


    vec3 t = vec3(0, 0, 0);
    vec3 s = vec3(1.0, 1.0, 1.0);
    vec3 r = vec3(0, 0, 0);
    

    std::vector<Color> m_materials;

    Color lightColor = Color(1, 1, 1);
    vec3 light = vec3(0, 0, 2.5);
    ImGuiIO io;
    float tick = 0;
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
