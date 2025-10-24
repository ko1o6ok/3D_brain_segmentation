// Необходимые инклюды
// -------------------
// ITK
// -------------------
// Работа с изображениями
#include <itkImage.h>
#include <itkImageFileReader.h>
// Конкретно для работы с .gipl файлами
#include <itkGiplImageIOFactory.h>

// Фильтры

// Для удаления шума
// #include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkMedianImageFilter.h>
// Для выделения границ
#include <itkGradientMagnitudeImageFilter.h>
// Для нормировки
#include <itkRescaleIntensityImageFilter.h>
// Для сегментации
#include <itkWatershedImageFilter.h>
// #include <itkConnectedThresholdImageFilter.h>
// #include <itkScalarImageKmeansImageFilter.h>
// #include <itkOtsuMultipleThresholdsImageFilter.h>
#include <itkLabelToRGBImageFilter.h>
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
#include "Parameters.h"


using InputPixelType = short;
using PixelType = float;
constexpr unsigned int Dimension = 3; // Исправлено: было Dimention
using InputImageType = itk::Image<InputPixelType, Dimension>;

using ImageType = itk::Image<PixelType, Dimension>;
using LabelImageType = itk::Image<unsigned long, Dimension>;
using RGBPixelType = itk::RGBPixel<unsigned char>;
using RGBImageType = itk::Image<RGBPixelType, Dimension>;
using LabelToRGBFilter = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;

int main()
{
    Params par;
    // Регистрируем поддержку GIPL формата
    itk::GiplImageIOFactory::RegisterOneFactory();
    
    std::cout << "=== ЗАПУСК VOLUME RENDERING ===" << std::endl;
    
    // ЭТАП 1: ЧТЕНИЕ ФАЙЛА
    std::cout << "1. Чтение файла..." << std::endl;
    
    using ReaderType = itk::ImageFileReader<InputImageType>;
    ReaderType::Pointer my_reader = ReaderType::New();
    my_reader->SetFileName("MRIcrop-orig.gipl");

    try {
        my_reader->Update();
        std::cout << "Файл прочитан успешно!" << std::endl;
        
        InputImageType::Pointer my_image = my_reader->GetOutput();
        InputImageType::RegionType region = my_image->GetLargestPossibleRegion();
        InputImageType::SizeType size = region.GetSize();
        std::cout << "   Размер изображения: " << size[0] << " x " << size[1] << " x " << size[2] << std::endl;
        
    } catch (itk::ExceptionObject &ex) {
        std::cerr << "Ошибка чтения файла: " << ex << std::endl;
        return 1;
    }

    // ЭТАП 2: ПАРАЛЛЕЛЬНАЯ ОБРАБОТКА - ОСНОВНОЕ ИЗОБРАЖЕНИЕ И СЕГМЕНТАЦИЯ
    std::cout << "2. Параллельная обработка..." << std::endl;

    // ПОДГОТОВКА ОСНОВНОГО ИЗОБРАЖЕНИЯ (ч/б)
    std::cout << "2.1. Подготовка основного изображения..." << std::endl;
    
    using ConnectorType = itk::ImageToVTKImageFilter<InputImageType>;
    auto connector = ConnectorType::New();
    connector->SetInput(my_reader->GetOutput());
    connector->Update();
    
    vtkImageData* volumeData = connector->GetOutput();
    
    // ДИАГНОСТИКА ОСНОВНЫХ ДАННЫХ
    int* dims = volumeData->GetDimensions();
    std::cout << "   Размеры основного объема: " << dims[0] << " x " << dims[1] << " x " << dims[2] << std::endl;
    
    double* range = volumeData->GetScalarRange();
    std::cout << "   Диапазон значений: " << range[0] << " - " << range[1] << std::endl;

    // ПОДГОТОВКА СЕГМЕНТАЦИИ
    std::cout << "2.2. Подготовка сегментации..." << std::endl;
    
    try {
        // ПАЙПЛАЙН СЕГМЕНТАЦИИ
        using CastFilterType = itk::CastImageFilter<InputImageType, ImageType>;
        auto castFilter = CastFilterType::New();
        castFilter->SetInput(my_reader->GetOutput());
        
        using MedianFilterType = itk::MedianImageFilter<ImageType, ImageType>;
        auto medianFilter = MedianFilterType::New();
        medianFilter->SetInput(castFilter->GetOutput());
        medianFilter->SetRadius(par.median_filter_medianRadius);
        
        using GradientFilterType = itk::GradientMagnitudeImageFilter<ImageType, ImageType>;
        auto gradientFilter = GradientFilterType::New();
        gradientFilter->SetInput(medianFilter->GetOutput());
        
        using WatershedFilterType = itk::WatershedImageFilter<ImageType>;
        auto watershedFilter = WatershedFilterType::New();
        watershedFilter->SetInput(gradientFilter->GetOutput());
        watershedFilter->SetThreshold(par.watershed_threshold);
        watershedFilter->SetLevel(par.watershed_level);
        
        using RGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
        auto rgbFilter = RGBFilterType::New();
        rgbFilter->SetInput(watershedFilter->GetOutput());
        rgbFilter->SetBackgroundValue(0);
        
        rgbFilter->Update();
        std::cout << " Watershed сегментация завершена" << std::endl;



        // Сегментация WATERSHED
        std::cout << "   Сегментация watershed..." << std::endl;

        // Проверим выход watershed (метки)
        LabelImageType::Pointer labelImage = watershedFilter->GetOutput();
        labelImage->Update();

        // Посчитаем уникальные метки
        using LabelIteratorType = itk::ImageRegionConstIterator<LabelImageType>;
        LabelIteratorType labelIt(labelImage, labelImage->GetLargestPossibleRegion());

        std::set<unsigned long> uniqueLabels;
        unsigned long minLabel = std::numeric_limits<unsigned long>::max();
        unsigned long maxLabel = 0;

        for (labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt) {
            unsigned long label = labelIt.Get();
            uniqueLabels.insert(label);
            if (label < minLabel) minLabel = label;
            if (label > maxLabel) maxLabel = label;
        }

        std::cout << "   - Уникальных меток: " << uniqueLabels.size() << std::endl;
        std::cout << "   - Диапазон меток: " << minLabel << " - " << maxLabel << std::endl;

        // Проверим выход RGB фильтра
        RGBImageType::Pointer rgbImage = rgbFilter->GetOutput();
        using RGBIteratorType = itk::ImageRegionConstIterator<RGBImageType>;
        RGBIteratorType rgbIt(rgbImage, rgbImage->GetLargestPossibleRegion());

        std::set<std::vector<unsigned char>> uniqueColors;
        for (rgbIt.GoToBegin(); !rgbIt.IsAtEnd(); ++rgbIt) {
            RGBPixelType pixel = rgbIt.Get();
            std::vector<unsigned char> color = {pixel[0], pixel[1], pixel[2]};
            uniqueColors.insert(color);
        }

        std::cout << "   - Уникальных цветов: " << uniqueColors.size() << std::endl;
        
        // КОННЕКТОР ДЛЯ СЕГМЕНТАЦИИ
        using SegmentationConnectorType = itk::ImageToVTKImageFilter<RGBImageType>;
        auto segmentation_connector = SegmentationConnectorType::New();
        segmentation_connector->SetInput(rgbFilter->GetOutput());
        segmentation_connector->Update();
        
        vtkImageData* segmentationData = segmentation_connector->GetOutput();

        // ПРЕОБРАЗУЕМ RGB В RGBA (добавляем альфа-канал)
        std::cout << "   Преобразование RGB в RGBA..." << std::endl;

        vtkSmartPointer<vtkImageData> rgbaData = vtkSmartPointer<vtkImageData>::New();
        rgbaData->SetDimensions(segmentationData->GetDimensions());
        rgbaData->AllocateScalars(VTK_UNSIGNED_CHAR, 4); // 4 компонента: RGBA

        int* dims = segmentationData->GetDimensions();
        for (int z = 0; z < dims[2]; z++) {
            for (int y = 0; y < dims[1]; y++) {
                for (int x = 0; x < dims[0]; x++) {
                    unsigned char* rgb = static_cast<unsigned char*>(
                        segmentationData->GetScalarPointer(x, y, z));
                    unsigned char* rgba = static_cast<unsigned char*>(
                        rgbaData->GetScalarPointer(x, y, z));
                    
                    // Копируем RGB
                    rgba[0] = rgb[0]; // R
                    rgba[1] = rgb[1]; // G  
                    rgba[2] = rgb[2]; // B
                    rgba[3] = 255;    // A - полностью непрозрачный
                }
            }
        }

        // Используем RGBA данные вместо RGB
        segmentationData = rgbaData;
        std::cout << " Данные преобразованы в RGBA" << std::endl;
        
        // ДИАГНОСТИКА СЕГМЕНТАЦИИ
        std::cout << "   Диагностика сегментации:" << std::endl;
        int* segDims = segmentationData->GetDimensions();
        std::cout << "   - Размеры: " << segDims[0] << " x " << segDims[1] << " x " << segDims[2] << std::endl;
        
        double* segRange = segmentationData->GetScalarRange();
        std::cout << "   - Диапазон скаляров: " << segRange[0] << " - " << segRange[1] << std::endl;
        
        int numComponents = segmentationData->GetNumberOfScalarComponents();
        std::cout << "   - Количество компонент: " << numComponents << std::endl;
        
        // ЭТАП 3: СОЗДАНИЕ VOLUME RENDERING ДЛЯ ОБОИХ ИЗОБРАЖЕНИЙ
        std::cout << "3. Создание Volume Rendering..." << std::endl;
        
        // 3.1. ОСНОВНОЙ VOLUME (ч/б мозг)
        std::cout << "3.1. Настройка основного volume..." << std::endl;
        
        vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volumeMapper = 
            vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
        volumeMapper->SetInputData(volumeData);
        volumeMapper->SetSampleDistance(0.3);
        volumeMapper->SetImageSampleDistance(1.5);
        
        vtkSmartPointer<vtkVolumeProperty> volumeProperty = 
            vtkSmartPointer<vtkVolumeProperty>::New();
        volumeProperty->ShadeOff();
        volumeProperty->SetInterpolationTypeToLinear();
        
        // Функция непрозрачности для основного изображения
        vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = 
            vtkSmartPointer<vtkPiecewiseFunction>::New();
        compositeOpacity->AddPoint(0,     0.00 * par.brainOpacity);
        compositeOpacity->AddPoint(20,    0.00 * par.brainOpacity);
        compositeOpacity->AddPoint(40,    0.15 * par.brainOpacity);
        compositeOpacity->AddPoint(70,    0.25 * par.brainOpacity);
        compositeOpacity->AddPoint(100,   0.40 * par.brainOpacity);
        compositeOpacity->AddPoint(150,   0.65 * par.brainOpacity);
        compositeOpacity->AddPoint(200,   0.80 * par.brainOpacity);
        compositeOpacity->AddPoint(269,   0.90 * par.brainOpacity);
        
        volumeProperty->SetScalarOpacity(compositeOpacity);
        
        // Цветовая функция - ч/б
        vtkSmartPointer<vtkColorTransferFunction> color = 
            vtkSmartPointer<vtkColorTransferFunction>::New();
        color->AddRGBPoint(0,     0.00, 0.00, 0.00);
        color->AddRGBPoint(40,    0.15, 0.15, 0.15);
        color->AddRGBPoint(70,    0.35, 0.35, 0.35);
        color->AddRGBPoint(100,   0.55, 0.55, 0.55);
        color->AddRGBPoint(150,   0.75, 0.75, 0.75);
        color->AddRGBPoint(200,   0.90, 0.90, 0.90);
        color->AddRGBPoint(269,   1.00, 1.00, 1.00);
        
        volumeProperty->SetColor(color);
        
        vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
        volume->SetMapper(volumeMapper);
        volume->SetProperty(volumeProperty);
        
        // 3.2. VOLUME ДЛЯ СЕГМЕНТАЦИИ с vtkGPUVolumeRayCastMapper
        std::cout << "3.2. Настройка volume для сегментации (GPU mapper)..." << std::endl;

        vtkSmartPointer<vtkGPUVolumeRayCastMapper> segmentationMapper = 
            vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
        segmentationMapper->SetInputData(segmentationData);
        segmentationMapper->SetAutoAdjustSampleDistances(0);
        segmentationMapper->SetSampleDistance(0.5);

        vtkSmartPointer<vtkVolumeProperty> segmentationProperty = 
            vtkSmartPointer<vtkVolumeProperty>::New();
        segmentationProperty->ShadeOff();
        segmentationProperty->SetInterpolationTypeToLinear();
        segmentationProperty->IndependentComponentsOff();  // Важно!

        // Функция непрозрачности
        vtkSmartPointer<vtkPiecewiseFunction> segmentationOpacity = 
            vtkSmartPointer<vtkPiecewiseFunction>::New();
        segmentationOpacity->AddPoint(0, 0.0 * par.segmentationOpacity);
        segmentationOpacity->AddPoint(1, 0.4 * par.segmentationOpacity);
        segmentationOpacity->AddPoint(255, 0.6 * par.segmentationOpacity);

        segmentationProperty->SetScalarOpacity(segmentationOpacity);

        vtkSmartPointer<vtkVolume> segmentationVolume = vtkSmartPointer<vtkVolume>::New();
        segmentationVolume->SetMapper(segmentationMapper);
        segmentationVolume->SetProperty(segmentationProperty);
        
        std::cout << "  Оба volume созданы" << std::endl;
        
        // ЭТАП 4: СОЗДАНИЕ СЦЕНЫ
        std::cout << "4. Создание сцены..." << std::endl;
        
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        
        // ПРАВИЛЬНЫЙ ПОРЯДОК: сначала основной volume, потом сегментация поверх
        renderer->AddVolume(volume);                    // Основной ч/б мозг
        renderer->AddVolume(segmentationVolume);        // Цветная сегментация поверх
        
        renderer->SetBackground(0.1, 0.1, 0.2);
        renderer->ResetCamera();
        
        // Настройка камеры
        vtkCamera* camera = renderer->GetActiveCamera();
        camera->Azimuth(30);
        camera->Elevation(30);
        
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
        
        interactor->Start();
        
    } catch (itk::ExceptionObject &ex) {
        std::cerr << "Ошибка в пайплайне сегментации: " << ex << std::endl;
        return 1;
    }

    std::cout << "Программа завершена" << std::endl;
    return 0;
}