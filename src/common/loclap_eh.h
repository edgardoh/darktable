#pragma once
/*
    This file is part of darktable,
    copyright (c) 2016 johannes hanika.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

void loclap_eh_internal(
    const float *const input,   // input buffer in some Labx or yuvx format
    float *const out,           // output buffer with colour
    const int wd,               // width and
    const int ht,               // height of the input buffer
    const float sigma,          // user param: separate shadows/midtones/highlights
    const float shadows,        // user param: lift shadows
    const float highlights,     // user param: compress highlights
    const float clarity,        // user param: increase clarity/local contrast
    const float denoise,        // user param: increase denoise
    const int use_sse2);        // switch on sse optimised version, if available

void loclap_eh(
    const float *const input,   // input buffer in some Labx or yuvx format
    float *const out,           // output buffer with colour
    const int wd,               // width and
    const int ht,               // height of the input buffer
    const float sigma,          // user param: separate shadows/midtones/highlights
    const float shadows,        // user param: lift shadows
    const float highlights,     // user param: compress highlights
    const float clarity,        // user param: increase clarity/local contrast
		const float denoise)        // user param: increase denoise
{
  loclap_eh_internal(input, out, wd, ht, sigma, shadows, highlights, clarity, denoise, 0);
}

#if defined(__SSE2__)
void loclap_eh_sse2(
    const float *const input,   // input buffer in some Labx or yuvx format
    float *const out,           // output buffer with colour
    const int wd,               // width and
    const int ht,               // height of the input buffer
    const float sigma,          // user param: separate shadows/midtones/highlights
    const float shadows,        // user param: lift shadows
    const float highlights,     // user param: compress highlights
    const float clarity,        // user param: increase clarity/local contrast
		const float denoise)        // user param: increase denoise
{
  loclap_eh_internal(input, out, wd, ht, sigma, shadows, highlights, clarity, denoise, 1);
}
#endif
