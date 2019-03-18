
#ifndef DT_DEVELOP_LOCALLAPLACIAN_EH_H
#define DT_DEVELOP_LOCALLAPLACIAN_EH_H

void loclap_LocalLaplacian(float *im_in, float *im_out, const int width, const int height, 
                            const float alpha, const float beta, const float *const scales, const float preview_scale, const int use_sse);

#endif

