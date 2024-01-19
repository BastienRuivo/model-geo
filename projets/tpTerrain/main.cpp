
//! \file tuto7_camera.cpp reprise de tuto7.cpp mais en derivant AppCamera, avec gestion automatique d'une camera.


#include "wavefront.h"
#include "uniforms.h"
#include "texture.h"

#include "orbiter.h"
#include "draw.h"        
#include "app_camera.h"        // classe Application a deriver

#include "HeightField.h"

#include "FastNoiseLite.h"

#include <random>

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


Point randomPointOnHemisphere(float u, float v) {
    float cos_theta = u;
    float sin_theta = sqrt(1 - u * u);
    float phi = 2 * M_PI * v;
    return RotationX(90)(Point(cos(phi) * sin_theta, cos_theta, sin_theta * sin(phi)));
}

class TP : public AppCamera
{
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
        //glewInit();

        ImGui::StyleColorsDark();

        noiser.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
        noiser.SetFractalType(FastNoiseLite::FractalType_FBm);
        noiser.SetFractalOctaves(2);
        noiser.SetFractalLacunarity(2.0);
        noiser.SetFractalGain(0.5);


        // decrire un repere / grille 
        tw = 1024;
        th = 1024;
        m_terrain = HeightField(Point(0, 0, 0), Point(10, 1.0, 10), tw, th);

        m_program = read_program("data/shaders/mesh.glsl");

        m_kirby = read_texture(0, "data/kirby.jpg");


        for (size_t i = 0; i < m_terrain.getWidth(); i++)
        {
            for (size_t j = 0; j < m_terrain.getHeight(); j++)
            {
                double x = (double)i / (double)m_terrain.getWidth();
                double y = (double)j / (double)m_terrain.getHeight();
                double value = noiser.GetNoise(x * tw, y * th);
                // clamp value between 0 and 1
                value = (value + 1) * 0.5;
                m_terrain(i, j) = value;
            }
        }

        //m_terrain.loadFromFile("data/terrain/heightmap.png", Point(0, 0, 0), Point(10, 2.5, 10));

        ImageData img = m_terrain.toImage();
        
        // generate texture
        glGenTextures(1, &m_texture);
        glBindTexture(GL_TEXTURE_2D, m_texture);
        // rgb texture
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img.width, img.height, 0, GL_RGB, GL_UNSIGNED_BYTE, img.pixels.data());
        glGenerateMipmap(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, 0);


        //m_terrain.smooth();
        m_terrain.polygonize(tw, th);

        m_grid= make_grid(10);
        

        Point pmin, pmax;
        m_terrain.mesh.bounds(pmin, pmax);
        pmin = RotationX(-90)(pmin);
        pmax = RotationX(-90)(pmax);
        camera().lookat(pmin, pmax);
        // std::cout<<"{";
        // for (size_t i = 0; i < 64; i++)
        // {
        //     auto p = randomPointOnHemisphere(distribution(generator), distribution(generator));

        //     // Mesh line = Mesh(GL_LINES);
        //     // line.color(1, 0, 0);
        //     // line.color(1, 0, 0);
        //     // line.vertex(Point(0, 0, 0));
        //     // line.vertex(p);
        //     // m_lines.push_back(line);

        //     std::cout<<"Point("<<p.x<<", "<<p.y<<", "<<p.z<<"), ";
        // }
        // std::cout<<"}"<<std::endl;
        // for(float i = 0; i < 10; i+=0.5) {
        //     for(float j = 0; j < 10; j+=0.5) {
        //         Point p = Point(i, m_terrain.Height(i, j), j);
        //         float sky = m_terrain.accesibility(p);
        //     }
        // }
        


        
        
        // etat openGL par defaut
        glClearColor(0.2f, 0.2f, 0.2f, 1.f);        // couleur par defaut de la fenetre
        
        glClearDepth(1.f);                          // profondeur par defaut
        glDepthFunc(GL_LESS);                       // ztest, conserver l'intersection la plus proche de la camera
        glEnable(GL_DEPTH_TEST);                    // activer le ztest

        return 0;   // pas d'erreur, sinon renvoyer -1
    }
    void updateImage() {
        // update texture
        glBindTexture(GL_TEXTURE_2D, m_texture);
        // rgb texture
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_img.width, m_img.height, 0, GL_RGB, GL_UNSIGNED_BYTE, m_img.pixels.data());
        glGenerateMipmap(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    // destruction des objets de l'applicatEEion
    int quit( )
    {
        //m_repere.release();
        return 0;   // pas d'erreur
    }

    void doUI() {
        ImGui_ImplSdlGL3_NewFrame(m_window);
        //ImGui::NewFrame();

        // render your GUI
        SDL_Event event;
        SDL_PollEvent(&event);
        
        ImGui_ImplSdlGL3_ProcessEvent(&event);
        ImGui::Begin("UI");


        if(ImGui::Checkbox("Wireframe", &isWireframe)) {
            if(isWireframe)
                glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
            else
                glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        }
        
        // start a frame
        if(ImGui::Button("Img To Terrain")) {
            m_terrain.UpdateFromImage(m_img);
            m_terrain.polygonize(256, 256);
        }
        ImGui::SameLine();
        if(ImGui::Button("Terrain to Img")) {
            m_img = m_terrain.toImage();
            updateImage();
        }
        ImGui::SameLine();
        if(ImGui::Button("Reset")) {
            m_terrain.loadFromFile("data/terrain/heightmap.png", Point(0, 0, 0), Point(10, 0.25, 10));
            m_terrain.polygonize(tw, th);
        }
        if(ImGui::Button("Slope")) {
            m_img = m_terrain.Slope(tw, th);
            updateImage();
        }
        ImGui::SameLine();
        if(ImGui::Button("Gradient Magnitude")) {
            m_img = m_terrain.GradientMagnitude(tw, th);
            updateImage();
        }
        ImGui::SameLine();
        if(ImGui::Button("Laplacian")) {
            m_img = m_terrain.Laplacian(tw, th);
            updateImage();
        }
        ImGui::SameLine();
        if(ImGui::Button("Accessible")) {
            m_img = m_terrain.Accessible(tw, th);
            updateImage();
        }
        ImGui::SameLine();
        if(ImGui::Button("Normal Color")) {
            m_img = m_terrain.NormalColor(tw, th);
            updateImage();
        }
        ImGui::Image((void*)(intptr_t)m_texture, ImVec2(512, 512));
        ImGui::End();

        ImGui::Render();
        ImGui_ImplSdlGL3_RenderDrawData(ImGui::GetDrawData());
    }
    
    void doKeyboard() {
        if (io.WantCaptureKeyboard) return;
    }

    // dessiner une nouvelle image
    int render( )
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        draw(m_grid, Identity(), camera());
        // Dessinne le repere et la lumiere

        Transform model= Identity();
        Transform view= camera().view();
        Transform projection= camera().projection(window_width(), window_height(), 45);

        glUseProgram(m_program);
        draw(m_terrain.mesh, Identity(), camera(), m_texture);

        for (size_t i = 0; i < m_terrain.lines.size(); i++)
        {
            draw(m_terrain.lines[i], Identity(), camera());
        }

        doUI();
        doKeyboard();
        return 1;
    }

protected:
    HeightField m_terrain;
    Mesh m_grid;
    GLuint m_texture;
    GLuint m_kirby;
    ImageData m_img;

    int tw, th;

    bool isWireframe = false;

    std::vector<Mesh> m_lines;

    GLuint m_program;

    ImGuiIO io;

    FastNoiseLite noiser;

    // default random engine
    std::default_random_engine generator;
    // random distribution
    std::uniform_real_distribution<float> distribution{0.0, 1.0};

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
