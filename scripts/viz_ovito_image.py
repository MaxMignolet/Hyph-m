from ovito.io import import_file
from ovito import modifiers
from ovito.vis import Viewport

def setup_particle_types(frame, data):
    types = data.particles_.particle_types_
    # types.type_by_id_(1).name = "Cr"
    types.type_by_id_(1).radius = 0.645
    # types.type_by_id_(2).name = "I"
    types.type_by_id_(2).radius = 0.695

for imode in range(20,24): # for CrI3
# for imode in range(15,17): # for CrBr3
  pipeline = import_file("outdata/dump_"+f'{imode:02}'+"_*")
  pipeline.modifiers.append(modifiers.WrapPeriodicImagesModifier())
  pipeline.add_to_scene()


  pipeline.modifiers.append(setup_particle_types)

  vp = Viewport()
  # vp.camera_pos = (0,0,10)
  pixel_size = (1920,1080) # 1080p
  # top view
  vp.camera_dir = (0,0,-1)
  vp.zoom_all()
  # vp.render_anim(filename="dir_viz/animation_"+f'{imode:02}'+"_top.mp4")
  for iframe in range(len(pipeline.frames)):
    vp.render_image(filename="dir_viz/image_"+f'{imode:02}'+'_'+f'{iframe:04}'+"_top.png",frame=iframe,size=pixel_size,crop=True,alpha=False)


  vp.camera_dir = (1,1,-0.4)
  vp.zoom_all()
  # vp.render_anim(filename="dir_viz/animation_"+f'{imode:02}'+"_front.mp4")
  for iframe in range(len(pipeline.frames)):
    vp.render_image(filename="dir_viz/image_"+f'{imode:02}'+'_'+f'{iframe:04}'+"_front.png",frame=iframe,size=pixel_size,crop=True,alpha=False)
  pipeline.remove_from_scene()
