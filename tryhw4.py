from uwimg import *

a = load_image("data/dog_a.jpg")
b = load_image("data/dog_b.jpg")
flow = optical_flow_images(b, a, 15, 8)
draw_flow(a, flow, 8)
save_image(a, "lines")

c = load_image("data/dog_a.jpg")
bf = box_filter_image(c, 8)
save_image(bf, "boxfilter")

#optical_flow_webcam(15,4,8)
