/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iomanip>
#include <iostream>
#include <cstdlib>

#include "mvs/settings.h"
#include "mvs/dmrecon.h"
#include "core/scene.h"
#include "core/view.h"
#include "util/timer.h"
#include "util/arguments.h"
#include "util/system.h"
#include "util/tokenizer.h"
#include "util/file_system.h"



struct AppSettings
{
    std::string scene_path;
    std::string ply_dest = "recon";
    int master_id = -1;
    std::vector<int> view_ids;
    int max_pixels = 1500000;
    bool force_recon = false;
    bool write_ply = false;
    mvs::Settings mvs;
};



int
main (int argc, char** argv)
{
    if(argc<3){
        std::cout<<"usage: scendir scale"<<std::endl;
        return -1;
    }

    AppSettings conf;

    // 场景文件夹
    conf.scene_path = argv[1];
    // 获取图像尺度
    std::stringstream stream1(argv[2]);
    stream1>>conf.mvs.scale;

    /* Load MVE scene. */
    core::Scene::Ptr scene;
    try
    {
        scene = core::Scene::create(conf.scene_path);
        scene->get_bundle();
    }
    catch (std::exception& e)
    {
        std::cerr << "Error loading scene: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    /* Settings for Multi-view stereo */
    conf.mvs.writePlyFile = conf.write_ply;
    conf.mvs.plyPath = util::fs::join_path(conf.scene_path, conf.ply_dest);

     core::Scene::ViewList& views(scene->get_views());
     if (conf.view_ids.empty())
     {
        std::cout << "Reconstructing all views..." << std::endl;
        for(std::size_t i = 0; i < views.size(); ++i)
             conf.view_ids.push_back(i);
     }
     else
     {
         std::cout << "Reconstructing views from list..." << std::endl;
     }

     // 对于每一个视角单独进行重建
     util::WallTimer timer;
     for (std::size_t i = 0; i < conf.view_ids.size(); ++i)
     {
            std::size_t id = conf.view_ids[i];
            if (id >= views.size())
            {
                std::cout << "Invalid ID " << id << ", skipping!" << std::endl;
                continue;
            }

            if (views[id] == nullptr || !views[id]->is_camera_valid())
                continue;

            /* Setup MVS. */
            mvs::Settings settings(conf.mvs);
            settings.refViewNr = id;

            std::string embedding_name = "depth-L" + util::string::get(settings.scale);
            if (!conf.force_recon && views[id]->has_image(embedding_name))
                continue;

            try {
                // 重建场景
                mvs::DMRecon recon(scene, settings);
                recon.start();
                views[id]->save_view();
            }
            catch (std::exception &err)
            {
                std::cerr << err.what() << std::endl;
            }
     }


    std::cout << "Reconstruction took "<< timer.get_elapsed() << "ms." << std::endl;

    /* Save scene */
    std::cout << "Saving views back to disc..." << std::endl;scene->save_views();

    return EXIT_SUCCESS;
}
