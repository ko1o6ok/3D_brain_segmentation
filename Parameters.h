// Parameters.h
#pragma once

struct Params {
    // Параметры watershed
    double watershed_threshold = 0.005;
    double watershed_level = 0.2;

    // Параметры фильтра шумов
    int median_filter_medianRadius = 2;

    // Прозрачность
    double brainOpacity = 1.0;
    double segmentationOpacity = 0.04;
};