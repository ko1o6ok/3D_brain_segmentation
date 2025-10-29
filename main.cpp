// Необходимые инклюды
// -------------------
// ITK
// -------------------
// Работа с изображениями
#include <itkImage.h>
#include <itkImageFileReader.h>
// Конкретно для работы с .gipl файлами
// #include <itkGiplImageIOFactory.h>
#include <itkMetaImageIOFactory.h>
#include <itkMetaImageIO.h>
// Фильтры

// Для удаления шума
// #include <itkGradientAnisotropicDiffusionImageFilter.h>
// #include <itkMedianImageFilter.h>
#include <itkCurvatureFlowImageFilter.h>
// Для выделения границ
#include <itkGradientMagnitudeImageFilter.h>
#include <itkCannyEdgeDetectionImageFilter.h>
// Для нормировки
#include <itkMinimumMaximumImageCalculator.h>
#include <itkRescaleIntensityImageFilter.h>
// Для сегментации
// #include <itkWatershedImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkComposeImageFilter.h>
#include <itkVectorConfidenceConnectedImageFilter.h>

#include <vtkImageReslice.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
// #include <itkScalarImageKmeansImageFilter.h>
// #include <itkOtsuMultipleThresholdsImageFilter.h>
// #include <itkLabelToRGBImageFilter.h>
// Переход к VTK
#include <itkImageToVTKImageFilter.h>


// -------------------
// VTK
// -------------------
// Умные указатели
#include <vtkSmartPointer.h>
// Для добавления акторов в сцену
#include <vtkImageActor.h>
// Для сцены
#include <vtkRenderer.h>
// Камера
#include <vtkCamera.h>
// Окно
#include <vtkRenderWindow.h>
// Интерактор для окна
#include <vtkRenderWindowInteractor.h>
// Взаимодействие типа TrackBall
#include <vtkInteractorStyleTrackballCamera.h>
// Для Volume Rendering - УПРОЩЁННАЯ ВЕРСИЯ
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkGPUVolumeRayCastMapper.h>


#include <vtkVolumeProperty.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolume.h>


#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>

// Рисование меток, где находятся сиды
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
// Для сообщений
#include <iostream>
#include <itkCastImageFilter.h>

// MATH
// #include <vector>
// #include <cmath>

// Рефакторинг
// #include "Parameters.h"


using InputPixelType = short;
using PixelType = float;
constexpr unsigned int Dimension = 3; // Исправлено: было Dimention
using InputImageType = itk::Image<InputPixelType, Dimension>;

using ImageType = itk::Image<PixelType, Dimension>;
using VectorImageType = itk::Image<itk::Vector<PixelType, 2>, Dimension>;
using LabelImageType = itk::Image<unsigned char, Dimension>;
using RGBPixelType = itk::RGBPixel<unsigned char>;
using RGBImageType = itk::Image<RGBPixelType, Dimension>;
// using LabelToRGBFilter = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;

// ДОБАВЬТЕ ЭТУ ФУНКЦИЮ ПЕРЕД main():
// Функция для создания визуальных маркеров сидов
void AddSeedMarkers(vtkRenderer* renderer, ImageType::Pointer image, 
                   const std::vector<std::array<int, 3>>& seeds,
                   const std::vector<std::array<double, 3>>& colors) {
    
    ImageType::SpacingType spacing = image->GetSpacing();
    ImageType::PointType origin = image->GetOrigin();
    
    for (size_t i = 0; i < seeds.size(); ++i) {
        // Преобразуем индекс в физические координаты
        ImageType::PointType physicalPoint;
        ImageType::IndexType idx = {{seeds[i][0], seeds[i][1], seeds[i][2]}};
        image->TransformIndexToPhysicalPoint(idx, physicalPoint);
        
        // Создаем сферу-маркер
        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetCenter(physicalPoint[0], physicalPoint[1], physicalPoint[2]);
        sphere->SetRadius(2.0); // Размер маркера (в мм)
        sphere->SetPhiResolution(12);
        sphere->SetThetaResolution(12);
        
        // Мэппер для сферы
        vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        sphereMapper->SetInputConnection(sphere->GetOutputPort());
        
        // Актор для сферы
        vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
        sphereActor->SetMapper(sphereMapper);
        sphereActor->GetProperty()->SetColor(colors[i][0], colors[i][1], colors[i][2]);
        sphereActor->GetProperty()->SetOpacity(0.8);
        
        // Добавляем в сцену
        renderer->AddActor(sphereActor);
        
        std::cout << "  Маркер " << i+1 << ": (" << seeds[i][0] << ", " << seeds[i][1] << ", " << seeds[i][2] 
                  << ") -> (" << physicalPoint[0] << ", " << physicalPoint[1] << ", " << physicalPoint[2] << ")" << std::endl;
    }
}


int main()
{
    // Параметры анизотропной диффузии
    // int num_iter = 5; // Кол-во итераций
    // float time_step = 0.03; // Шаг по времени
    // float conductance = 3.0; // Чувствительность к границам

    // Параметры watershed
    // float threshold = 0.005; //0.005// минимальный порог, ниже которого пиксели
    // // не рассматриваются как потенциальные маркеры или области для сегментации
    // float level = 0.2; //0.2 // Уровень заполнения

    // Параметры медианного фильтра
    // int medianRadius = 2;

    // Регистрируем поддержку GIPL формата
    // itk::GiplImageIOFactory::RegisterOneFactory();
    itk::MetaImageIOFactory::RegisterOneFactory();
    
    std::cout << "=== ЗАПУСК VOLUME RENDERING ===" << std::endl;
    
    // ЭТАП 1: ЧТЕНИЕ ФАЙЛА
    std::cout << "1. Чтение файла..." << std::endl;
    
    using ReaderType = itk::ImageFileReader<InputImageType>;
    ReaderType::Pointer my_reader = ReaderType::New();
    my_reader->SetFileName("patient_109_mr_T2.mhd");

    ReaderType::Pointer my_reader_add = ReaderType::New();
    my_reader_add->SetFileName("patient_109_mr_T1.mhd");
    // itk::MetaImageIO::Pointer metaIO = itk::MetaImageIO::New();
    // metaIO->SetByteOrderToBigEndian();
    // my_reader->SetImageIO(metaIO);

    // Загрузка файла с T2
    try {
        my_reader->Update();
        std::cout << "Файл .mhd прочитан успешно!" << std::endl;
        InputImageType::Pointer my_image = my_reader->GetOutput();
        InputImageType::RegionType region = my_image->GetLargestPossibleRegion();
        InputImageType::SizeType size = region.GetSize();
        InputImageType::SpacingType spacing = my_image->GetSpacing();
        std::cout << "   Размер изображения: " << size[0] << " x " << size[1] << " x " << size[2] << std::endl;
        std::cout << "   Пространственное разрешение: " << spacing[0] << " x " << spacing[1] << " x " << spacing[2] << std::endl;
        // ДИАГНОСТИКА: проверяем реальный диапазон данных
        itk::ImageRegionConstIterator<InputImageType> it(my_image, region);
        it.GoToBegin();
        InputPixelType minVal = itk::NumericTraits<InputPixelType>::max();
        InputPixelType maxVal = itk::NumericTraits<InputPixelType>::min();
        while (!it.IsAtEnd()) {
            InputPixelType val = it.Get();
            if (val < minVal) minVal = val;
            if (val > maxVal) maxVal = val;
            ++it;
        }
        std::cout << "   Реальный диапазон данных: " << minVal << " - " << maxVal << std::endl;
    } catch (itk::ExceptionObject &ex) {
        std::cerr << "Ошибка чтения .mhd файла: " << ex << std::endl;
        return 1;
    }

    // Загрузка файла с T1
    try {
        my_reader_add->Update();
        std::cout << "Файл .mhd прочитан успешно!" << std::endl;
        InputImageType::Pointer my_image_add = my_reader_add->GetOutput();
        InputImageType::RegionType region_add = my_image_add->GetLargestPossibleRegion();
        InputImageType::SizeType size_add = region_add.GetSize();
        InputImageType::SpacingType spacing = my_image_add->GetSpacing();
        std::cout << "   Размер изображения: " << size_add[0] << " x " << size_add[1] << " x " << size_add[2] << std::endl;
        std::cout << "   Пространственное разрешение: " << spacing[0] << " x " << spacing[1] << " x " << spacing[2] << std::endl;
        // ДИАГНОСТИКА: проверяем реальный диапазон данных
        itk::ImageRegionConstIterator<InputImageType> it(my_image_add, region_add);
        it.GoToBegin();
        InputPixelType minVal = itk::NumericTraits<InputPixelType>::max();
        InputPixelType maxVal = itk::NumericTraits<InputPixelType>::min();
        while (!it.IsAtEnd()) {
            InputPixelType val = it.Get();
            if (val < minVal) minVal = val;
            if (val > maxVal) maxVal = val;
            ++it;
        }
        std::cout << "   Реальный диапазон данных: " << minVal << " - " << maxVal << std::endl;
    } catch (itk::ExceptionObject &ex) {
        std::cerr << "Ошибка чтения .mhd файла: " << ex << std::endl;
        return 1;
    }

    // ЭТАП 2: ПАРАЛЛЕЛЬНАЯ ОБРАБОТКА - ОСНОВНОЕ ИЗОБРАЖЕНИЕ И СЕГМЕНТАЦИЯ
    std::cout << "2. Обработка..." << std::endl;
    
    // ПОДГОТОВКА СЕГМЕНТАЦИИ
    std::cout << "2.1. Подготовка сегментации..." << std::endl;
    
    try {
        // T2 - пайплайн
        // ПАЙПЛАЙН СЕГМЕНТАЦИИ
        using CastFilterType = itk::CastImageFilter<InputImageType, ImageType>;
        auto castFilter = CastFilterType::New();
        castFilter->SetInput(my_reader->GetOutput());

        using DenoiseFilterType = itk::CurvatureFlowImageFilter<ImageType, ImageType>;
        auto denoiseFilter = DenoiseFilterType::New();
        denoiseFilter->SetInput(castFilter->GetOutput());
        denoiseFilter->SetTimeStep(0.125);
        denoiseFilter->SetNumberOfIterations(5);


        using RescaleFilterType = itk::RescaleIntensityImageFilter<ImageType, ImageType>;
        auto rescaleFilter = RescaleFilterType::New();
        rescaleFilter->SetInput(denoiseFilter->GetOutput());
        rescaleFilter->SetOutputMinimum(0.0f);
        rescaleFilter->SetOutputMaximum(255.0f);
        rescaleFilter->Update();



        // T1 - пайплайн
        using CastFilterType = itk::CastImageFilter<InputImageType, ImageType>;
        auto castFilter_add = CastFilterType::New();
        castFilter_add->SetInput(my_reader_add->GetOutput());

        using DenoiseFilterType = itk::CurvatureFlowImageFilter<ImageType, ImageType>;
        auto denoiseFilter_add = DenoiseFilterType::New();
        denoiseFilter_add->SetInput(castFilter_add->GetOutput());
        denoiseFilter_add->SetTimeStep(0.125);
        denoiseFilter_add->SetNumberOfIterations(5);


        using RescaleFilterType = itk::RescaleIntensityImageFilter<ImageType, ImageType>;
        auto rescaleFilter_add = RescaleFilterType::New();
        rescaleFilter_add->SetInput(denoiseFilter_add->GetOutput());
        rescaleFilter_add->SetOutputMinimum(0.0f);
        rescaleFilter_add->SetOutputMaximum(255.0f);
        rescaleFilter_add->Update();

        // Композиция и сегментация серого вещества
        using ComposeFilterType = itk::ComposeImageFilter<ImageType, VectorImageType>;
        auto compose_filter = ComposeFilterType::New();
        compose_filter->SetInput1(rescaleFilter->GetOutput());
        compose_filter->SetInput2(rescaleFilter_add->GetOutput());

        int seed1_x = 165, seed1_y = 178, seed1_z = 26;
        int seed2_x = 98, seed2_y = 165, seed2_z = 26; 
        int seed3_x = 205, seed3_y = 125, seed3_z = 26;
        int  seed4_x = 173, seed4_y = 205, seed4_z = 26;

        using VectorConfidenceConnectedFilterType = itk::VectorConfidenceConnectedImageFilter<VectorImageType, LabelImageType>;
        auto conf_filter = VectorConfidenceConnectedFilterType::New();
        conf_filter->SetInput(compose_filter->GetOutput());
        conf_filter->SetNumberOfIterations(1);
        conf_filter->SetMultiplier(0.1);
        conf_filter->SetReplaceValue(255);
        conf_filter->AddSeed({seed1_x, seed1_y, seed1_z});
        conf_filter->AddSeed({seed2_x, seed2_y, seed2_z});
        conf_filter->AddSeed({seed3_x, seed3_y, seed3_z});
        conf_filter->AddSeed({seed4_x, seed4_y, seed4_z});

        std::cout << "   Добавлены сиды для VectorConfidenceConnected:" << std::endl;
        std::cout << "   Сид 1: (" << seed1_x << ", " << seed1_y << ", " << seed1_z << ")" << std::endl;
        std::cout << "   Сид 2: (" << seed2_x << ", " << seed2_y << ", " << seed2_z << ")" << std::endl;
        std::cout << "   Сид 3: (" << seed3_x << ", " << seed3_y << ", " << seed3_z << ")" << std::endl;
        std::cout << "   Сид 4: (" << seed4_x << ", " << seed4_y << ", " << seed4_z << ")" << std::endl;
        

        // Белое вещество сегментация
        int white_seed1_x = 150, white_seed1_y = 150, white_seed1_z = 26;

        int white_seed2_x = 100, white_seed2_y = 150, white_seed2_z = 26;

        int white_seed3_x = 152, white_seed3_y = 65, white_seed3_z = 26;

        int white_seed4_x = 120, white_seed4_y = 80, white_seed4_z = 26;

        std::cout << "2.3. Сегментация белого вещества (мультимодальная)..." << std::endl;

        using WhiteMatterFilterType = itk::VectorConfidenceConnectedImageFilter<VectorImageType, LabelImageType>;
        auto white_matter_filter = WhiteMatterFilterType::New();
        white_matter_filter->SetInput(compose_filter->GetOutput());
        white_matter_filter->SetNumberOfIterations(3);  // Больше итераций для лучшего роста
        white_matter_filter->SetMultiplier(2.0);        // Более строгий множитель для белого вещества
        white_matter_filter->SetReplaceValue(180);      // Уникальное значение для белого вещества
        // white_matter_filter->SetInitialNeighborhoodRadius(2);  // Больший радиус для лучшей статистики

        // Добавьте сиды
        white_matter_filter->AddSeed({white_seed1_x, white_seed1_y, white_seed1_z});
        white_matter_filter->AddSeed({white_seed2_x, white_seed2_y, white_seed2_z});
        white_matter_filter->AddSeed({white_seed3_x, white_seed3_y, white_seed3_z});
        white_matter_filter->AddSeed({white_seed4_x, white_seed4_y, white_seed4_z});

        std::cout << "   Добавлены сиды для белого вещества:" << std::endl;
        std::cout << "   Сид 1: (" << white_seed1_x << ", " << white_seed1_y << ", " << white_seed1_z << ")" << std::endl;
        std::cout << "   Сид 2: (" << white_seed2_x << ", " << white_seed2_y << ", " << white_seed2_z << ")" << std::endl;
        std::cout << "   Сид 3: (" << white_seed3_x << ", " << white_seed3_y << ", " << white_seed3_z << ")" << std::endl;
        std::cout << "   Сид 4: (" << white_seed4_x << ", " << white_seed4_y << ", " << white_seed4_z << ")" << std::endl;

        white_matter_filter->Update();



        // ДИАГНОСТИКА: проверяем, есть ли сегментированные пиксели
        LabelImageType::Pointer whiteResult = white_matter_filter->GetOutput();
        itk::ImageRegionConstIterator<LabelImageType> whiteIt(whiteResult, whiteResult->GetLargestPossibleRegion());
        whiteIt.GoToBegin();
        unsigned long whitePixelCount = 0;
        while (!whiteIt.IsAtEnd()) {
            if (whiteIt.Get() == 180) {
                whitePixelCount++;
            }
            ++whiteIt;
        }
        std::cout << "   Сегментировано пикселей белого вещества: " << whitePixelCount << std::endl;

        // СЕГМЕНТАЦИЯ ГЛАЗ - МИНИМАЛЬНЫЙ КОД
        std::cout << "2.2. Сегментация глаз..." << std::endl;

        using ConnectedFilterType = itk::ConnectedThresholdImageFilter<ImageType, ImageType>;
        auto eyeFilter = ConnectedFilterType::New();
        eyeFilter->SetInput(rescaleFilter->GetOutput());

        // ПРИМЕРНЫЕ КООРДИНАТЫ ГЛАЗ
        int right_eye_x = 100;
        int right_eye_y = 55;
        int right_eye_z = 20;
        eyeFilter->AddSeed({right_eye_x, right_eye_y, right_eye_z});

        int left_eye_x = 172;
        int left_eye_y = 55;
        int left_eye_z = 20; 
        eyeFilter->AddSeed({left_eye_x, left_eye_y, left_eye_z});

        // Диапазон интенсивностей для глаз
        eyeFilter->SetLower(80);   // 60 минимальная интенсивность глаз
        eyeFilter->SetUpper(200);  // 100 максимальная интенсивность глаз
        eyeFilter->SetReplaceValue(255); // Белый цвет для сегментированных глаз

        eyeFilter->Update();

        // Подключаем сегментированные глаза к VTK
        using EyeConnectorType = itk::ImageToVTKImageFilter<ImageType>;
        auto eye_connector = EyeConnectorType::New();
        eye_connector->SetInput(eyeFilter->GetOutput());
        eye_connector->Update();

        vtkImageData* eyeData = eye_connector->GetOutput();


        // Подключаем результат VectorConfidenceConnected к VTK
        using ConfidenceConnectorType = itk::ImageToVTKImageFilter<LabelImageType>;
        auto confidence_connector = ConfidenceConnectorType::New();
        confidence_connector->SetInput(conf_filter->GetOutput());
        confidence_connector->Update();

        vtkImageData* confidenceData = confidence_connector->GetOutput();

        using WhiteMatterConnectorType = itk::ImageToVTKImageFilter<LabelImageType>;
        auto white_matter_connector = WhiteMatterConnectorType::New();
        white_matter_connector->SetInput(white_matter_filter->GetOutput());
        white_matter_connector->Update();

        vtkImageData* whiteMatterData = white_matter_connector->GetOutput();






        // using SegmentationConnectorType = itk::ImageToVTKImageFilter<LabelImageType>;
        using SegmentationConnectorType = itk::ImageToVTKImageFilter<ImageType>;
        auto segmentation_connector = SegmentationConnectorType::New();
        segmentation_connector->SetInput(rescaleFilter->GetOutput());
        segmentation_connector->Update();

        vtkImageData* segmentationData = segmentation_connector->GetOutput();

        // Диагностика данных
        double* dataRange = segmentationData->GetScalarRange();
        int* dims = segmentationData->GetDimensions();
    
        std::cout << "   Диапазон интенсивностей после обработки: " << dataRange[0] << " - " << dataRange[1] << std::endl;
        std::cout << "   Размеры VTK данных: " << dims[0] << " x " << dims[1] << " x " << dims[2] << std::endl;
        
        // ЭТАП 3: СОЗДАНИЕ VOLUME RENDERING ДЛЯ ОБОИХ ИЗОБРАЖЕНИЙ
        std::cout << "3. Создание Volume Rendering..." << std::endl;
        
        // 3.2. VOLUME ДЛЯ СЕГМЕНТАЦИИ с vtkGPUVolumeRayCastMapper
        std::cout << "3.2. Настройка volume для сегментации (GPU mapper)..." << std::endl;

        vtkSmartPointer<vtkGPUVolumeRayCastMapper> segmentationMapper = 
            vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
        segmentationMapper->SetInputData(segmentationData);
        segmentationMapper->SetAutoAdjustSampleDistances(0);
        segmentationMapper->SetSampleDistance(0.3);

        vtkSmartPointer<vtkVolumeProperty> segmentationProperty = 
            vtkSmartPointer<vtkVolumeProperty>::New();
        segmentationProperty->ShadeOff();
        segmentationProperty->SetInterpolationTypeToLinear();
        segmentationProperty->IndependentComponentsOn();  // Важно!

        // Функция непрозрачности
        vtkSmartPointer<vtkPiecewiseFunction> segmentationOpacity = 
            vtkSmartPointer<vtkPiecewiseFunction>::New();

        segmentationOpacity->AddPoint(0,   0.0);
        segmentationOpacity->AddPoint(55,  0.0 ); //50 0.05
        
        segmentationOpacity->AddPoint(65,  1.0 * 0.03);
        segmentationOpacity->AddPoint(255, 1.0 * 0.03);

        vtkSmartPointer<vtkColorTransferFunction> segmentationColor = 
        vtkSmartPointer<vtkColorTransferFunction>::New();
    
        segmentationColor->AddRGBPoint(0,   0.0, 0.0, 0.0);   // Черный
        segmentationColor->AddRGBPoint(50,  0.5, 0.5, 0.5); // 0.2  // Темно-серый
        // segmentationColor->AddRGBPoint(120, 0.6, 0.6, 0.6);   // Серый
        // segmentationColor->AddRGBPoint(180, 0.7, 0.7, 0.7);   // Светло-серый
        segmentationColor->AddRGBPoint(80, 1.0, 1.0, 1.0); 
        segmentationColor->AddRGBPoint(255, 1.0, 1.0, 1.0);  // Белый

        segmentationProperty->SetScalarOpacity(segmentationOpacity);
        segmentationProperty->SetColor(segmentationColor);

        vtkSmartPointer<vtkVolume> segmentationVolume = vtkSmartPointer<vtkVolume>::New();
        segmentationVolume->SetMapper(segmentationMapper);
        segmentationVolume->SetProperty(segmentationProperty);
        
        std::cout << "  Volume создан" << std::endl;


        // VOLUME ДЛЯ СЕГМЕНТИРОВАННЫХ ГЛАЗ
        std::cout << "3.3. Настройка volume для глаз..." << std::endl;

        vtkSmartPointer<vtkGPUVolumeRayCastMapper> eyeMapper = 
            vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
        eyeMapper->SetInputData(eyeData);
        eyeMapper->SetSampleDistance(0.3);

        vtkSmartPointer<vtkVolumeProperty> eyeProperty = 
            vtkSmartPointer<vtkVolumeProperty>::New();
        eyeProperty->ShadeOff();
        eyeProperty->SetInterpolationTypeToLinear();

        // Глаза - КРАСНЫЕ и полупрозрачные
        vtkSmartPointer<vtkPiecewiseFunction> eyeOpacity = 
            vtkSmartPointer<vtkPiecewiseFunction>::New();
        eyeOpacity->AddPoint(0, 0.0);    // Фон - прозрачный
        eyeOpacity->AddPoint(255, 0.8);  //0.8 Глаза - полупрозрачные

        vtkSmartPointer<vtkColorTransferFunction> eyeColor = 
            vtkSmartPointer<vtkColorTransferFunction>::New();
        eyeColor->AddRGBPoint(0, 0.0, 0.0, 0.0);    // Фон - черный
        eyeColor->AddRGBPoint(255, 1.0, 0.0, 0.0);  // Глаза - КРАСНЫЕ

        eyeProperty->SetScalarOpacity(eyeOpacity);
        eyeProperty->SetColor(eyeColor);

        vtkSmartPointer<vtkVolume> eyeVolume = vtkSmartPointer<vtkVolume>::New();
        eyeVolume->SetMapper(eyeMapper);
        eyeVolume->SetProperty(eyeProperty);

        std::cout << "  Volume для глаз создан" << std::endl;




        // VOLUME ДЛЯ РЕЗУЛЬТАТА VectorConfidenceConnected
        std::cout << "3.4. Настройка volume для VectorConfidenceConnected..." << std::endl;

        vtkSmartPointer<vtkGPUVolumeRayCastMapper> confidenceMapper = 
            vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
        confidenceMapper->SetInputData(confidenceData);
        confidenceMapper->SetSampleDistance(0.3);

        vtkSmartPointer<vtkVolumeProperty> confidenceProperty = 
            vtkSmartPointer<vtkVolumeProperty>::New();
        confidenceProperty->ShadeOff();
        confidenceProperty->SetInterpolationTypeToLinear();

        // Результат VectorConfidenceConnected - СИНИЙ
        vtkSmartPointer<vtkPiecewiseFunction> confidenceOpacity = 
            vtkSmartPointer<vtkPiecewiseFunction>::New();
        confidenceOpacity->AddPoint(0, 0.0);    // Фон - прозрачный
        confidenceOpacity->AddPoint(255, 0.7);  // Сегментация - полупрозрачная

        vtkSmartPointer<vtkColorTransferFunction> confidenceColor = 
            vtkSmartPointer<vtkColorTransferFunction>::New();
        confidenceColor->AddRGBPoint(0, 0.0, 0.0, 0.0);    // Фон - черный
        confidenceColor->AddRGBPoint(255, 0.0, 1.0, 0.0);  // Сегментация - 

        confidenceProperty->SetScalarOpacity(confidenceOpacity);
        confidenceProperty->SetColor(confidenceColor);

        vtkSmartPointer<vtkVolume> confidenceVolume = vtkSmartPointer<vtkVolume>::New();
        confidenceVolume->SetMapper(confidenceMapper);
        confidenceVolume->SetProperty(confidenceProperty);

        std::cout << "  Volume для VectorConfidenceConnected создан" << std::endl;
        


        std::cout << "3.5. Настройка volume для белого вещества..." << std::endl;

        vtkSmartPointer<vtkGPUVolumeRayCastMapper> whiteMatterMapper = 
            vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
        whiteMatterMapper->SetInputData(whiteMatterData);
        whiteMatterMapper->SetSampleDistance(0.3);

        vtkSmartPointer<vtkVolumeProperty> whiteMatterProperty = 
            vtkSmartPointer<vtkVolumeProperty>::New();
        whiteMatterProperty->ShadeOff();
        whiteMatterProperty->SetInterpolationTypeToLinear();

        // Белое вещество - СИНИЙ цвет
        vtkSmartPointer<vtkPiecewiseFunction> whiteMatterOpacity = 
            vtkSmartPointer<vtkPiecewiseFunction>::New();
        whiteMatterOpacity->AddPoint(0, 0.0);    // Фон - прозрачный
        whiteMatterOpacity->AddPoint(180, 0.6);  // Белое вещество - полупрозрачное

        vtkSmartPointer<vtkColorTransferFunction> whiteMatterColor = 
            vtkSmartPointer<vtkColorTransferFunction>::New();
        whiteMatterColor->AddRGBPoint(0, 0.0, 0.0, 0.0);      // Фон - черный
        whiteMatterColor->AddRGBPoint(180, 0.0, 0.4, 1.0);    // Белое вещество - СИНИЙ

        whiteMatterProperty->SetScalarOpacity(whiteMatterOpacity);
        whiteMatterProperty->SetColor(whiteMatterColor);

        vtkSmartPointer<vtkVolume> whiteMatterVolume = vtkSmartPointer<vtkVolume>::New();
        whiteMatterVolume->SetMapper(whiteMatterMapper);
        whiteMatterVolume->SetProperty(whiteMatterProperty);

        std::cout << "  Volume для белого вещества создан" << std::endl;
        // ЭТАП 4: СОЗДАНИЕ СЦЕНЫ
        std::cout << "4. Создание сцены..." << std::endl;
        
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        
        // ПРАВИЛЬНЫЙ ПОРЯДОК: сначала основной volume, потом сегментация поверх                    // Основной ч/б мозг
        renderer->AddVolume(segmentationVolume);        // Цветная сегментация поверх
        
        renderer->AddVolume(confidenceVolume);

        renderer->AddVolume(whiteMatterVolume);

        renderer->AddVolume(eyeVolume);

        // ВИЗУАЛИЗАЦИЯ СИДОВ В 3D ОКНЕ
        std::cout << "4.1. Добавление маркеров сидов в 3D окно..." << std::endl;

        // Сиды для глаз
        // std::vector<std::array<int, 3>> eyeSeeds = {
        //     {right_eye_x, right_eye_y, right_eye_z},
        //     {left_eye_x, left_eye_y, left_eye_z}
        // };

        // std::vector<std::array<double, 3>> eyeColors = {
        //     {1.0, 0.0, 0.0},  // Красный - правый глаз
        //     {1.0, 0.5, 0.0}   // Оранжевый - левый глаз
        // };

        // AddSeedMarkers(renderer, tempImage, eyeSeeds, eyeColors);

        // Сиды для VectorConfidenceConnected
        std::vector<std::array<int, 3>> confSeeds = {
            {seed1_x, seed1_y, seed1_z},
            {seed2_x, seed2_y, seed2_z},
            {seed3_x, seed3_y, seed3_z},
            {seed4_x, seed4_y, seed4_z}
        };

        std::vector<std::array<double, 3>> confColors = {
            {0.0, 1.0, 0.0},  // Зеленый
            {0.0, 1.0, 1.0},  // Голубой  
            {1.0, 0.0, 1.0},   // Фиолетовый
            {1.0, 1.0, 0.0}
        };

        AddSeedMarkers(renderer, rescaleFilter->GetOutput(), confSeeds, confColors);

        std::vector<std::array<int, 3>> whiteMatterSeeds = {
            {white_seed1_x, white_seed1_y, white_seed1_z},
            {white_seed2_x, white_seed2_y, white_seed2_z},
            {white_seed3_x, white_seed3_y, white_seed3_z},
            {white_seed4_x, white_seed4_y, white_seed4_z}
        };

        std::vector<std::array<double, 3>> whiteMatterColors = {
            {0.0, 0.4, 1.0},  // Синий
            {0.2, 0.5, 1.0},  // Светло-синий
            {0.0, 0.3, 0.8},  // Темно-синий
            {0.4, 0.6, 1.0}   // Очень светлый синий
        };

        AddSeedMarkers(renderer, rescaleFilter->GetOutput(), whiteMatterSeeds, whiteMatterColors);



        
        renderer->SetBackground(0.0, 0.0, 0.0);
        renderer->ResetCamera();
        
        // Настройка камеры
        vtkCamera* camera = renderer->GetActiveCamera();
        camera->Azimuth(0);
        camera->Elevation(-90); // -90
        
        vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
        window->AddRenderer(renderer);
        window->SetSize(1000, 1000);
        window->SetWindowName("3D Volume Rendering - Мозг с сегментацией");
        window->SetMultiSamples(0);
        
        vtkSmartPointer<vtkRenderWindowInteractor> interactor = 
            vtkSmartPointer<vtkRenderWindowInteractor>::New();
        interactor->SetRenderWindow(window);
        
        vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
            vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
        interactor->SetInteractorStyle(style);
        
        std::cout << "5. Запуск рендеринга..." << std::endl;
        window->Render();
        
        std::cout << "Рендеринг запущен!" << std::endl;
        std::cout << "Управление:" << std::endl;
        std::cout << "  - ЛКМ: вращение" << std::endl;
        std::cout << "  - ПКМ: масштабирование" << std::endl;
        std::cout << "  - СКМ: перемещение" << std::endl;








        // СОЗДАНИЕ ОТДЕЛЬНОГО ОКНА С 2D СРЕЗОМ
        std::cout << "6. Создание окна с 2D срезом..." << std::endl;

        // Создаем аксиальный срез
        vtkSmartPointer<vtkImageReslice> sliceReslice = vtkSmartPointer<vtkImageReslice>::New();
        sliceReslice->SetInputData(segmentationData);
        sliceReslice->SetOutputDimensionality(2);
        sliceReslice->SetResliceAxesDirectionCosines(1,0,0, 0,1,0, 0,0,1);
        sliceReslice->SetResliceAxesOrigin(0, 0, 78); // Средний срез
        sliceReslice->Update();

        // Создаем актор для среза
        vtkSmartPointer<vtkImageActor> sliceActor = vtkSmartPointer<vtkImageActor>::New();
        sliceActor->GetMapper()->SetInputConnection(sliceReslice->GetOutputPort());

        // Создаем отдельный рендерер и окно для 2D среза
        vtkSmartPointer<vtkRenderer> sliceRenderer = vtkSmartPointer<vtkRenderer>::New();
        sliceRenderer->AddActor(sliceActor);
        sliceRenderer->SetBackground(0.2, 0.2, 0.2); // Темный фон для лучшего контраста


        // СОЗДАНИЕ СРЕЗА ДЛЯ СЕГМЕНТАЦИИ CONFIDENCE CONNECTED
        std::cout << "6.2. Создание среза для сегментации Confidence Connected..." << std::endl;

        // Создаем срез из результата VectorConfidenceConnected на том же уровне
        vtkSmartPointer<vtkImageReslice> confidenceSliceReslice = vtkSmartPointer<vtkImageReslice>::New();
        confidenceSliceReslice->SetInputData(confidenceData);
        confidenceSliceReslice->SetOutputDimensionality(2);
        confidenceSliceReslice->SetResliceAxesDirectionCosines(1,0,0, 0,1,0, 0,0,1);
        confidenceSliceReslice->SetResliceAxesOrigin(0, 0, 26 * 3); // Тот же уровень, что и у сидов
        confidenceSliceReslice->Update();

        // Создаем актор для среза сегментации
        vtkSmartPointer<vtkImageActor> confidenceSliceActor = vtkSmartPointer<vtkImageActor>::New();
        confidenceSliceActor->GetMapper()->SetInputConnection(confidenceSliceReslice->GetOutputPort());

        // Настраиваем внешний вид сегментации в 2D
        // Сегментация будет отображаться как ЗЕЛЕНЫЕ области поверх исходного изображения
        // confidenceSliceActor->SetOpacity(0.3); // Полупрозрачность для видимости исходного изображения под ней
        // confidenceSliceActor->GetProperty()->SetColor(0.0, 1.0, 0.0);
        // Создаем lookup table для зелёного цвета
        vtkSmartPointer<vtkLookupTable> confidenceLUT = vtkSmartPointer<vtkLookupTable>::New();
        confidenceLUT->SetNumberOfColors(256);
        confidenceLUT->SetTableRange(0, 255);
        confidenceLUT->Build();

        // Устанавливаем черный цвет для фона (значение 0)
        confidenceLUT->SetTableValue(0, 0.0, 0.0, 0.0, 0.0); // полностью прозрачный

        // Устанавливаем зелёный цвет для сегментированных областей (значение 255)
        confidenceLUT->SetTableValue(255, 0.0, 1.0, 0.0, 0.7); // зелёный, полупрозрачный

        // Применяем LUT к данным
        vtkSmartPointer<vtkImageMapToColors> confidenceColorMapper = vtkSmartPointer<vtkImageMapToColors>::New();
        confidenceColorMapper->SetLookupTable(confidenceLUT);
        confidenceColorMapper->SetInputConnection(confidenceSliceReslice->GetOutputPort());
        confidenceColorMapper->Update();

        // Обновляем актор с новыми данными
        confidenceSliceActor->GetMapper()->SetInputConnection(confidenceColorMapper->GetOutputPort());


        // Добавляем в сцену ПОСЛЕ исходного среза, но ДО маркеров сидов
        sliceRenderer->AddActor(confidenceSliceActor);

        std::cout << "Срез сегментации Confidence Connected добавлен (зеленый, полупрозрачный)" << std::endl;

        std::cout << "6.4. Создание среза для белого вещества..." << std::endl;

        vtkSmartPointer<vtkImageReslice> whiteMatterSliceReslice = vtkSmartPointer<vtkImageReslice>::New();
        whiteMatterSliceReslice->SetInputData(whiteMatterData);
        whiteMatterSliceReslice->SetOutputDimensionality(2);
        whiteMatterSliceReslice->SetResliceAxesDirectionCosines(1,0,0, 0,1,0, 0,0,1);
        whiteMatterSliceReslice->SetResliceAxesOrigin(0, 0, 26 * 3);
        whiteMatterSliceReslice->Update();

        // СОЗДАЕМ LUT ДЛЯ СИНЕГО ЦВЕТА (белое вещество)
        vtkSmartPointer<vtkLookupTable> whiteMatterLUT = vtkSmartPointer<vtkLookupTable>::New();
        whiteMatterLUT->SetNumberOfColors(256);
        whiteMatterLUT->SetTableRange(0, 255);
        whiteMatterLUT->Build();

        // Настраиваем цвета: 0 = прозрачный, 180 = синий
        for (int i = 0; i < 256; i++) {
            if (i == 180) {
                whiteMatterLUT->SetTableValue(i, 0.0, 0.4, 1.0, 0.6); // Синий, полупрозрачный
            } else {
                whiteMatterLUT->SetTableValue(i, 0.0, 0.0, 0.0, 0.0); // Полностью прозрачный
            }
        }

        // Применяем LUT к данным
        vtkSmartPointer<vtkImageMapToColors> whiteMatterColorMapper = vtkSmartPointer<vtkImageMapToColors>::New();
        whiteMatterColorMapper->SetLookupTable(whiteMatterLUT);
        whiteMatterColorMapper->SetInputConnection(whiteMatterSliceReslice->GetOutputPort());
        whiteMatterColorMapper->Update();

        // Создаем актор для среза белого вещества
        vtkSmartPointer<vtkImageActor> whiteMatterSliceActor = vtkSmartPointer<vtkImageActor>::New();
        whiteMatterSliceActor->GetMapper()->SetInputConnection(whiteMatterColorMapper->GetOutputPort());

        sliceRenderer->AddActor(whiteMatterSliceActor);

        // ДОБАВЛЯЕМ СИДЫ В 2D ОКНО
        std::cout << "6.1. Добавление маркеров сидов в 2D окно..." << std::endl;

        // Фильтруем сиды - оставляем только те, которые на текущем срезе
        std::vector<std::array<int, 3>> seedsOnSlice;
        std::vector<std::array<double, 3>> colorsOnSlice;

        // Добавляем сиды VectorConfidenceConnected (если они на этом срезе)  
        for (size_t i = 0; i < confSeeds.size(); ++i) {
                seedsOnSlice.push_back(confSeeds[i]);
                colorsOnSlice.push_back(confColors[i]);
        }

        // ДОБАВЛЯЕМ СИДЫ БЕЛОГО ВЕЩЕСТВА
        for (size_t i = 0; i < whiteMatterSeeds.size(); ++i) {
            seedsOnSlice.push_back(whiteMatterSeeds[i]);
            colorsOnSlice.push_back(whiteMatterColors[i]);
        }

        // Создаем маркеры для 2D окна
        // AddSeedMarkers(sliceRenderer, rescaleFilter->GetOutput(), seedsOnSlice, colorsOnSlice);

        // СОЗДАЕМ МАРКЕРЫ ДЛЯ 2D ОКНА С Z=0
        for (size_t i = 0; i < seedsOnSlice.size(); ++i) {
            // Преобразуем индекс в физические координаты
            ImageType::PointType physicalPoint;
            ImageType::IndexType idx = {{seedsOnSlice[i][0], seedsOnSlice[i][1], seedsOnSlice[i][2]}};
            rescaleFilter->GetOutput()->TransformIndexToPhysicalPoint(idx, physicalPoint);
            
            // Ключевое исправление: устанавливаем Z=0 для 2D плоскости
            physicalPoint[2] = 0.0;
            
            // Создаем сферу-маркер
            vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
            sphere->SetCenter(physicalPoint[0], physicalPoint[1], physicalPoint[2]);
            sphere->SetRadius(2.0); // Размер маркера (в мм)
            sphere->SetPhiResolution(12);
            sphere->SetThetaResolution(12);
            
            // Мэппер для сферы
            vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            sphereMapper->SetInputConnection(sphere->GetOutputPort());
            
            // Актор для сферы
            vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
            sphereActor->SetMapper(sphereMapper);
            sphereActor->GetProperty()->SetColor(colorsOnSlice[i][0], colorsOnSlice[i][1], colorsOnSlice[i][2]);
            sphereActor->GetProperty()->SetOpacity(0.8);
            
            // Добавляем в сцену
            sliceRenderer->AddActor(sphereActor);
            
            std::cout << "  Маркер " << i+1 << ": (" << seedsOnSlice[i][0] << ", " << seedsOnSlice[i][1] << ", " << seedsOnSlice[i][2] 
                    << ") -> (" << physicalPoint[0] << ", " << physicalPoint[1] << ", " << physicalPoint[2] << ")" << std::endl;
        }




        vtkSmartPointer<vtkRenderWindow> sliceWindow = vtkSmartPointer<vtkRenderWindow>::New();
        sliceWindow->AddRenderer(sliceRenderer);
        sliceWindow->SetSize(600, 600);
        sliceWindow->SetWindowName("2D Срез - Аксиальный вид");
        sliceWindow->SetPosition(1050, 100); // Позиция справа от основного окна

        // Создаем интерактор для 2D окна (ОБЩИЙ с основным окном!)
        vtkSmartPointer<vtkRenderWindowInteractor> sliceInteractor = 
            vtkSmartPointer<vtkRenderWindowInteractor>::New();
        sliceInteractor->SetRenderWindow(sliceWindow);

        // Запускаем рендеринг 2D окна
        sliceWindow->Render();

        std::cout << "   Окно с 2D срезом создано" << std::endl;
        
        interactor->Start();
        
    } catch (itk::ExceptionObject &ex) {
        std::cerr << "Ошибка в пайплайне сегментации: " << ex << std::endl;
        return 1;
    }

    std::cout << "Программа завершена" << std::endl;
    return 0;
}