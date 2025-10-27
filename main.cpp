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
using LabelImageType = itk::Image<unsigned long, Dimension>;
using RGBPixelType = itk::RGBPixel<unsigned char>;
using RGBImageType = itk::Image<RGBPixelType, Dimension>;
// using LabelToRGBFilter = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;

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
        
        // ЭТАП 4: СОЗДАНИЕ СЦЕНЫ
        std::cout << "4. Создание сцены..." << std::endl;
        
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        
        // ПРАВИЛЬНЫЙ ПОРЯДОК: сначала основной volume, потом сегментация поверх                    // Основной ч/б мозг
        renderer->AddVolume(segmentationVolume);        // Цветная сегментация поверх
        renderer->AddVolume(eyeVolume);

        std::cout << "Добавление визуальных маркеров для семян..." << std::endl;

        // Получаем информацию о геометрии изображения для преобразования координат
        ImageType::Pointer tempImage = rescaleFilter->GetOutput();
        ImageType::SpacingType spacing = tempImage->GetSpacing();
        ImageType::PointType origin = tempImage->GetOrigin();

        // Координаты семян (те же, что используются в сегментации)
        std::vector<ImageType::IndexType> seeds = {
            {right_eye_x, right_eye_y, right_eye_z},  
            {left_eye_x, left_eye_y, left_eye_z} 
        };

        // Цвета для маркеров (разные цвета для разных семян)
        // std::vector<std::array<double, 3>> colors = {
        //     {1.0, 0.0, 0.0},  // Красный - первое семя
        //     {0.0, 1.0, 0.0}   // Зеленый - второе семя
        // };

        // for (size_t i = 0; i < seeds.size(); ++i) {
        //     // Преобразуем индекс в физические координаты
        //     ImageType::PointType physicalPoint;
        //     tempImage->TransformIndexToPhysicalPoint(seeds[i], physicalPoint);
            
        //     // Создаем сферу-маркер
        //     vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        //     sphere->SetCenter(physicalPoint[0], physicalPoint[1], physicalPoint[2]);
        //     sphere->SetRadius(2.0); // Размер маркера (в мм)
        //     sphere->SetPhiResolution(20);
        //     sphere->SetThetaResolution(20);
            
        //     // Мэппер для сферы
        //     vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        //     sphereMapper->SetInputConnection(sphere->GetOutputPort());
            
        //     // Актор для сферы
        //     vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
        //     sphereActor->SetMapper(sphereMapper);
        //     sphereActor->GetProperty()->SetColor(colors[i][0], colors[i][1], colors[i][2]); // Цвет
        //     sphereActor->GetProperty()->SetOpacity(0.8); // Полупрозрачность
            
        //     // Добавляем в сцену
        //     renderer->AddActor(sphereActor);
            
        //     std::cout << "  Маркер " << i+1 << " в координатах: (" 
        //             << seeds[i][0] << ", " << seeds[i][1] << ", " << seeds[i][2] << ")" 
        //             << " -> физические: (" << physicalPoint[0] << ", " 
        //             << physicalPoint[1] << ", " << physicalPoint[2] << ")" << std::endl;
        // }
        
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
        sliceReslice->SetResliceAxesOrigin(0, 0, 80); // Средний срез
        sliceReslice->Update();

        // Создаем актор для среза
        vtkSmartPointer<vtkImageActor> sliceActor = vtkSmartPointer<vtkImageActor>::New();
        sliceActor->GetMapper()->SetInputConnection(sliceReslice->GetOutputPort());

        // Создаем отдельный рендерер и окно для 2D среза
        vtkSmartPointer<vtkRenderer> sliceRenderer = vtkSmartPointer<vtkRenderer>::New();
        sliceRenderer->AddActor(sliceActor);
        sliceRenderer->SetBackground(0.2, 0.2, 0.2); // Темный фон для лучшего контраста

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