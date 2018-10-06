clear all; close all; clc
u0_color = imread('1.bmp');
u0_color = double(u0_color);
u0 = u0_color(:,:,1);
phi = ActiveContours_const(u0, 0.6 * 255 * 255,1, 1,1000);