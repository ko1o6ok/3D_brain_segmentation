# 3D Brain Segmentation using ITK + VTK

## Автор
Щербаков Павел, группа 3825М1ПМвм1

![Результаты сегментации](SegmentationExample.png)

## Описание
Проект реализует алгоритмы 3D-сегментации данных компьютерной томограммы для выделения следующих структур:
- Глаза
- Серое вещество
- Белое вещество

## Функциональность
- Обработка 3D данных томографии
- Сегментация на основе алгоритмов ITK
- Визуализация результатов с использованием VTK

Реализована возможность включения/выключения частей сегментации по нажатию на кнопки ON/OFF.
Во вспомогательном окне размещена информация о сегментированном срезе с наложенным исходным срезом.

## Технологии
- **Язык программирования**: C++
- **Библиотеки**:
  - [ITK](https://itk.org/) — для обработки и сегментации медицинских изображений
  - [VTK](https://vtk.org/) — для визуализации 3D-данных

## Структура проекта
Весь проект реализован в одном файле main.cpp

## Весь пайплайн обработки
```mermaid
graph TD
    subgraph "1. Загрузка данных (ITK)"
        A1["ImageFileReader: T2 .mhd"] --> A2["CastImageFilter: short → float"]
        A2 --> A3["CurvatureFlowImageFilter: удаление шума"]
        A3 --> A4["RescaleIntensityImageFilter: 0..255"]

        B1["ImageFileReader: T1 .mhd"] --> B2["CastImageFilter: short → float"]
        B2 --> B3["CurvatureFlowImageFilter: удаление шума"]
        B3 --> B4["RescaleIntensityImageFilter: 0..255"]
    end

    subgraph "2. Композиция и сегментация (ITK)"
        A4 --> C["ComposeImageFilter: T2 + T1 → VectorImage"]
        B4 --> C

        C --> D["VectorConfidenceConnectedImageFilter: серое вещество"]
        C --> E["VectorConfidenceConnectedImageFilter: белое вещество"]
        A4 --> F["ConnectedThresholdImageFilter: глаза"]
    end

    subgraph "3. Переход в VTK (ITK → VTK)"
        D --> D1["ImageToVTKImageFilter → vtkImageData: серое вещество"]
        E --> E1["ImageToVTKImageFilter → vtkImageData: белое вещество"]
        F --> F1["ImageToVTKImageFilter → vtkImageData: глаза"]
    end

    subgraph "4. Volume Rendering (VTK)"
        D1 --> D2["vtkGPUVolumeRayCastMapper: серое вещество"]
        E1 --> E2["vtkGPUVolumeRayCastMapper: белое вещество"]
        F1 --> F2["vtkGPUVolumeRayCastMapper: глаза"]

        D2 --> D3["vtkVolume: серое вещество"]
        E2 --> E3["vtkVolume: белое вещество"]
        F2 --> F3["vtkVolume: глаза"]
    end

    subgraph "5. Surface Rendering (VTK)"
        D1 --> D4["vtkContourFilter: изоповерхность 128"]
        E1 --> E4["vtkContourFilter: изоповерхность 90"]
        F1 --> F4["vtkContourFilter: изоповерхность 128"]

        D4 --> D5["vtkWindowedSincPolyDataFilter: сглаживание"]
        E4 --> E5["vtkWindowedSincPolyDataFilter: сглаживание"]
        F4 --> F5["vtkWindowedSincPolyDataFilter: сглаживание"]

        D5 --> D6["vtkPolyDataMapper"]
        E5 --> E6["vtkPolyDataMapper"]
        F5 --> F6["vtkPolyDataMapper"]

        D6 --> D7["vtkActor: контур серого"]
        E6 --> E7["vtkActor: контур белого"]
        F6 --> F7["vtkActor: контур глаз"]
    end

    subgraph "6. Визуализация и интерактивность (VTK)"
        D3 --> G{VTK Renderer}
        E3 --> G
        F3 --> G

        D7 --> G
        E7 --> G
        F7 --> G

        G --> H["VTK Render Window: 3D View"]
        G --> I["VTK Render Window: 2D Slice View"]

        J["AddSeedMarkers: сиды в 3D и 2D"] --> G
        K["MouseCallback: переключение видимости"] --> H
    end

    subgraph "7. 2D Срезы (VTK)"
        D1 --> L["vtkImageReslice: 2D срез серого"]
        E1 --> M["vtkImageReslice: 2D срез белого"]
        F1 --> N["vtkImageReslice: 2D срез глаз"]

        L --> L1["vtkImageActor: 2D срез"]
        M --> M1["vtkImageActor: 2D срез (цвет)"]
        N --> N1["vtkImageActor: 2D срез"]

        L1 --> I
        M1 --> I
        N1 --> I
    end

    style A1 fill:#e1f5fe
    style B1 fill:#e1f5fe
    style D fill:#f3e5f5
    style E fill:#f3e5f5
    style F fill:#f3e5f5
    style D7 fill:#a5d6a7
    style E7 fill:#90caf9
    style F7 fill:#ef9a9a
    style D3 fill:#a5d6a7
    style E3 fill:#90caf9
    style F3 fill:#ef9a9a
```
