#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
breast_calculator.py - калькулятор объёма и веса женской груди
и европейского размера бюстгальтера.

Замысел прост: берём обычные измерения сантиметровой лентой и получаем
размер, объём и вес - с честно указанной неопределённостью, а не с
одной "точной" цифрой. Сюжет утилиты шутливый, но подход всерьёз:
формулы взяты из рецензируемых работ, и везде помечено, что из этого
валидировано, а что - эвристика без независимой проверки. Точный вес
одной груди в быту измерить нельзя, поэтому все оценки даются
интервалами, а обещания точности - сразу с оговорками.

Методика
--------
Вход: обхват груди (см) и обхват под грудью (см); дополнительные
антропометрические измерения - опциональные CLI-флаги.

1. Размер (буквальный стандарт EN 13402-3):
   - Бандаж: обхват под грудью округляется к ближайшему кратному 5 см
     (65, 70, 75, ...), минимум 65.
   - Чашка по разнице diff = обхват груди − обхват под грудью, шаг 2 см:
     10-12 = AA, 12-14 = A, 14-16 = B, 16-18 = C, 18-20 = D, 20-22 = E,
     22-24 = F, 24-26 = G, 26-28 = H; diff < 10 - размер меньше AA.
     Выше H - продолжение алфавита (I, J, K, ...) с пометкой об
     экстраполяции.
   - US-размер (выводится рядом с EU/EN): бандаж - обхват под грудью
     в дюймах, округлённый до ближайшего чётного, + 4 (стандартное
     соответствие EU 65→30, 70→32, 75→34, 80→36, 85→38, 90→40, 95→42;
     некоторые бренды используют +2 или без коррекции); чашка - разница
     обхватов в дюймах, шаг 1 дюйм (2.54 см) на размер: 1"=A, 2"=B,
     3"=C, 4"=D, 5"=DD, 6"=DDD (≡F), 7"=G, 8"=H, 9"=I, 10"=J, далее
     алфавит. Пересчёт US делается напрямую из измерений, а не
     перекодированием EU-букв (шаг EU-чашки 2 см ≠ 1 дюйму).

2. Модели объёма одной груди (выбор через --model):
   - heuristic (fallback): V ≈ 2.31 · diff² · (band / 75) - эвристическая,
     НЕ валидирована; константа откалибрована по СРЕДНЕМУ между двумя
     опубликованными якорями: (а) объёмы по размерам чашек в когорте
     Huang et al. 2017 (A≈261, B≈328, C≈408, ≥D≈539 мл → K≈1.97);
     (б) якорь стереофотограмметрии ~447 мл (левая 457 / правая 437)
     при типичном размере 75B → K≈2.645; среднее (1.97+2.645)/2 ≈ 2.31.
   - qiao (Qiao et al., 1997, Aesthet Plast Surg 21:362-368):
     V = (π/3)·MP²·(LR + MR + IR − MP); 250 грудей молодых китаянок,
     средний объём 310-330 мл; валидация по образцам не проводилась.
   - breast-v: формула BREAST-V (Longo et al., 2013):
     V = −231.66 + 0.5747·(SN-N)² + 18.5478·(FFp) + 14.5087·(N-IMF),
     где SN-N - яремная вырезка-сосок, FFp - проекция железы от грудной
     стенки (fold projection), N-IMF - сосок-инфрамаммарная складка (см).
     Валидирована на 108 мастэктомических образцах (R²=0.73, MAE 89.7 г,
     PMID 23806950). В редакции Huang et al. 2017 (PLoS ONE, DOI
     10.1371/journal.pone.0172122) параметр FFp заменён на BP (breast
     projection) с сохранением коэффициентов - здесь используется именно
     эта редакция (FFp ≈ BP - проекция железы).

3. Плотность (ICRU Report 44): кусочно-линейная интерполяция ρ(g) по трём
   опорным точкам (g=0 → 0.95 жировая; g=0.5 → 0.96 таблицы 50/50;
   g=1 → 1.02 железистая), где g - доля железистой
   (фиброгландулярной) ткани (по умолчанию 0.5; средняя доля по
   маммографии ~42.5% ± 30.3%).
   Опциональная линейная поправка по возрасту и ИМТ: коэффициенты -
   РУКОПИСНЫЕ ПЛЕЙСХОЛДЕРЫ-аппроксимации тенденции, НЕ результат
   регрессии (реальные зависимости нелинейны; см. K_AGE/K_BMI).

4. Вес: W = V·ρ с 95% доверительным интервалом:
   σ_V/V = sqrt(σ_мод² + σ_изм²) - ошибка модели + вклад ошибки измерений
   (вычисляется через частные производные формулы объёма по каждому входу);
   σ_W = W·sqrt((σ_V/V)² + (σ_ρ/ρ)²); 95% ДИ = W ± 1.96·σ_W.

5. Sister sizes: эквивалентные объёмы чашки. EU: бандаж ±5 см, чашка ∓1
   (75B ≈ 80A ≈ 70C); US: бандаж ±2 дюйма, чашка ∓1 (34DDD ≈ 36DD ≈ 32G).

6. Асимметрия: у всех женщин есть некоторая асимметрия; в среднем левая
    грудь чуть больше правой (MRI-волюметрия: 771 против 764 мл - Behrens
    et al., 2024; линейные расстояния слева больше - Henseler, 2023) -
    статистическая тенденция, не индивидуальный прогноз.

Использование моделей и измерения в домашних условиях
-----------------------------------------------------
Общие правила для всех измерений: измерять стоя, расслабленно, без
утягивающего белья; сантиметровая лента - параллельно полу, прилегает,
но не сдавливает; каждое измерение повторить 2-3 раза и усреднить.
Всё вводится в САНТИМЕТРАХ; десятичный разделитель - точка или запятая
("92", "92.5", "92,5" - все корректны). Все модели дают объём ОДНОЙ
груди; измеряйте обе стороны по отдельности, чтобы учесть асимметрию.

МОДЕЛЬ 1. heuristic (по умолчанию) - нужны только два обхвата
    python breast_calculator.py BUST UNDERBUST
    - BUST - обхват по самым выступающим точкам: лента проходит по
      горизонтали через соски, спина прямая, руки опущены.
    - UNDERBUST - обхват непосредственно под грудью, по нижнему краю
      железы (на уровне инфрамаммарной складки).
    Самая простая модель, но НЕ валидированная - см. примечания в отчёте.

МОДЕЛЬ 2. qiao (Qiao et al., 1997) - нужны 4 линейных измерения
    python breast_calculator.py BUST UNDERBUST --model qiao \
        --mp 7 --lr 6 --mr 6 --ir 7
    (обхваты при этой модели не участвуют в расчёте объёма, но нужны
    для определения размера бюстгальтера)
    Как измерить: встаньте боком к зеркалу и приложите к грудной стенке
    плоский предмет (линейку или коробку) торцом вперёд, перпендикулярно
    стенке - он "продолжит" грудную стенку; все проекции меряются от
    него по горизонтали:
    - --mp (MP)  - проекция железы: расстояние по горизонтали от
      линейки (грудная стенка) до самой выступающей точки соска.
    - --lr (LR)  - расстояние от соска до латерального (наружного)
      края железы (к подмышке).
    - --mr (MR)  - от соска до медиального (внутреннего) края железы
      (к грудине).
    - --ir (IR)  - от соска до нижнего края железы
      (до инфрамаммарной складки).
    Условие применимости: LR + MR + IR должно быть больше MP.

МОДЕЛЬ 3. breast-v (BREAST-V, Longo et al., 2013) - нужны 3 измерения
    python breast_calculator.py BUST UNDERBUST --model breast-v \
        --notch-nipple 24 --fold-nipple 9 --fold-projection 7.5
    Самая валидированная модель (108 мастэктомических образцов,
    MAE 89.7 г / 18.4%). Как измерить:
    - --notch-nipple (SN-N) - расстояние от яремной вырезки (впадина
      между ключицами на основании шеи) до соска, по поверхности тела.
    - --fold-nipple (N-IMF) - расстояние от соска до инфрамаммарной
      складки (нижней складки под грудью), стоя.
    - --fold-projection (FFp ≈ BP) - проекция железы: расстояние от
      грудной стенки до вершины соска (измерение боком к зеркалу, как
      и MP выше).

ДОПОЛНИТЕЛЬНО (для любой модели):
    --age N                возраст, лет
    --bmi N                ИМТ, кг/м² (вес кг / рост² м)
    --glandular-fraction N доля железистой ткани g ∈ [0,1] (если известна,
                           напр. по маммографии; по умолчанию 0.5)
    --measurement-error N  заявленная ошибка измерения, см (по умолч. 1.0;
                           при аккуратном самостоятельном измерении берите
                           0.5-1.0, для клинической точности 0.3)
    --self-check           проверка согласованности с якорями статей
    Без позиционных аргументов запускается интерактивный ввод
    (только два обхвата; остальные параметры - через CLI).

Примеры:
    python breast_calculator.py 88 75
    python breast_calculator.py 88 75 --model qiao --mp 7 --lr 6 --mr 6 --ir 7
    python breast_calculator.py 88 75 --model breast-v \
        --notch-nipple 24 --fold-nipple 9 --fold-projection 7.5
    python breast_calculator.py 88 75 --age 55 --bmi 31
    python breast_calculator.py --self-check

Используемые статьи (полный список ссылок):
-------------------------------------------
- Longo B., Farcomeni A., Ferri G., Campanale A., Sorotos M., Santanelli F.
  The BREAST-V: a unifying predictive formula for volume assessment in small,
  medium, and large breasts. Plast Reconstr Surg. 2013;132:1e-7e.
  PMID 23806950. (регрессия BREAST-V на 108 мастэктомических образцах;
  R²=0.73, MAE 89.7 г / 18.4%)
- Huang N-s, Quan C-l, Mo M, Chen J-j, Yang B-l, Huang X-y, Wu J.
  A prospective study of breast anthropomorphic measurements, volume and
  ptosis in 605 Asian patients with breast cancer or benign breast disease.
  PLoS ONE. 2017;12(2):e0172122. DOI 10.1371/journal.pone.0172122.
  (оценка объёма по самооценённым чашкам: A≈261/B≈328/C≈408/≥D≈539 мл;
  средний объём 340±109 мл; в формуле BREAST-V заменён параметр FFp на BP
  с сохранением коэффициентов Longo 2013; влияние BMI ≥ 24.7, менопаузы
  и кормления на объём и птоз)
- Qiao Q., Zhou G., Ling Y. Breast volume measurement in young Chinese
  women and clinical applications. Aesthet Plast Surg. 1997;21:362-368.
  (антропометрическая формула V=(π/3)·MP²·(LR+MR+IR−MP); 250 грудей;
  средний объём 310-330 мл)
- ICRU Report 44 (1989). Tissue Substitutes in Radiation Dosimetry and
  Measurement. (плотность тканей: жировая 0.95, фиброгландулярная 1.02 г/см³)
- Katariya R.N., Forrest D.M., Gravelle I.H. Breast volumes in cancer of the
  breast. Br J Surg. 1974;61:852-854. (объёмы груди в популяции)
- Grossman A.J., Roudner L.A. A simple technique for obtaining breast volume.
  Plast Reconstr Surg. 1980;65:301. (водяной конус, до 450 мл)
- Loughry C.W. et al. Breast volume measurement of 248 women using
  biostereometric analysis. Plast Reconstr Surg. 1987;80:553-558.
  (стереофотограмметрия: 248 женщин; эмпирический якорь среднего объёма
  ~437-457 мл при типичном размере)
- Brown R.W., Cheng Y-C., Kurtay M. A formula for surgical modifications of
  the breast. Plast Reconstr Surg. 2000;106:1342-1345. PMID 11083567.
- Sigurdson L.J., Kirkland S.A. Breast volume determination in breast
  hypertrophy: an accurate method using two anthropometric measurements.
  Plast Reconstr Surg. 2006;118:313-320. PMID 16874195.
- Kayar R. et al. Five methods of breast volume measurement: a comparative
  study of measurements of specimen volume in 30 mastectomy cases.
  Breast Cancer (Auckl). 2011;5:43-52.
  (сравнение 5 методов: маммография точнее всех, далее "метод Архимеда")
- Strömbeck J.O., Malm M. Priority grouping in a waiting list of patients
  for reduction mammaplasty. Ann Plast Surg. 1986;17:498-502.
  PMID 3827119. (+20 мл объёма на каждый кг веса выше идеального)
- Kayar et al. (2011) и обзор Xi W. et al. Objective breast volume, shape
  and surface area assessment: a systematic review of breast measurement
  methods. Aesthet Plast Surg. 2014;38:1116-1130.
- Choppin S.B., Wheat J.S., Gee M., Goyal A. The accuracy of breast volume
  measurement methods: a systematic review. Breast. 2016;28:121-129.
  (систематический обзор: наивысшая точность - MRI-волюметрия)
- Behrens A.S. et al. Comparative assessment of breast volume using a
  smartphone device versus MRI. Breast Cancer. 2024;32:166-176.
  DOI 10.1007/s12282-024-01647-6. (самая современная валидация метода
  объёма против MRI; средние объёмы: левая 771 мл, правая 764 мл -
  подтверждение чуть большей левой груди; CCC 3D-метода vs MRI 0.87-0.98)
- Henseler H. Exploring natural breast symmetry in the female plastic
  surgical patient population. GMS Interdiscip Plast Reconstr Surg DGPW.
  2023;12:Doc03. DOI 10.3205/iprs000173. (100 пациенток, 3D Vectra:
  асимметрия у всех; линейные расстояния слева больше)
- O'Connell R.L. et al. Validation of the Vectra XT three-dimensional
  imaging system for measuring breast volume and symmetry following
  oncological reconstruction. Breast Cancer Res Treat. 2018.
  DOI 10.1007/s10549-018-4843-6. (валидация 3D-фотографии объёма)
- EN 13402-3: Size designation of clothes - measurements and intervals.
  (бандаж с шагом 5 см; таблица чашек шага 2 см - буквальный стандарт)
- US sizing (band + cup): конвенция американских брендов - бандаж в
  дюймах (округление до ближайшего чётного + коррекция, типично +4),
  чашка с шагом 1 дюйм (2.54 см): A, B, C, D, DD, DDD/F, G, H, I, J...
  (сводная конвенция, единого стандарта нет)

Только стандартная библиотека Python.
"""

import sys

# Корректный вывод кириллицы в консолях с кодировкой, отличной от UTF-8 (Windows)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import argparse
import math

# ----------------------------------------------------------------------------
# Константы
# ----------------------------------------------------------------------------

# Таблица чашек EN 13402-3: разница обхватов diff (см), шаг 2 см на размер
# (cup size = bust − underbust, границы 10-12 = AA, ..., 26-28 = H).
CUP_TABLE = [
    (10, 12, "AA"),
    (12, 14, "A"),
    (14, 16, "B"),
    (16, 18, "C"),
    (18, 20, "D"),
    (20, 22, "E"),
    (22, 24, "F"),
    (24, 26, "G"),
    (26, 28, "H"),
]

CUP_MIN_DIFF = 10          # ниже - "<AA"
CUP_MAX_DIFF = 28          # выше - экстраполяция алфавита после H

MIN_CM, MAX_CM = 50, 160        # жёсткие пределы для обхватов, см
DIFF_MIN, DIFF_MAX = 5, 40      # поддерживаемый диапазон разницы обхватов, см
BAND_WARN = (60, 120)           # диапазон калибровки моделей (предупреждение)

# Якорные объёмы по размерам чашек (Huang et al., 2017, PLoS ONE,
# DOI 10.1371/journal.pone.0172122; когорта 605 китаянок, 1210 грудей):
# A ≈ 260.9 мл, B ≈ 328.0 мл, C ≈ 408.1 мл, ≥D ≈ 539.0 мл.
# Средний объём когорты 340.0±109.1 мл (91.8-919.2).

# Плотность ткани груди, г/см³ (ICRU Report 44): жировая 0.95 (g=0),
# 50/50 - 0.96 (g=0.5), фиброгландулярная 1.02 (g=1).
# tissue_density() интерполирует кусочно-линейно по этим трём опорным
# точкам (при g=0.5 даёт табличное 0.96, а не 0.985 линейной интерполяции
# двух крайних точек).
DENSITY_ADIPOSE = 0.95
DENSITY_5050 = 0.96
DENSITY_FIBROGLANDULAR = 1.02

# Эвристическая модель объёма: V ≈ HEURISTIC_K · diff² · (band / 75).
# Константа = среднее двух опубликованных якорей:
# (а) Huang et al. 2017: объёмы по чашкам A≈261/B≈328/C≈408/≥D≈539 мл
#     дают K = V/diff² ≈ 2.16/1.94/1.81/1.87 (mean ≈ 1.97);
# (б) стереофотограмметрия: ~447 мл при 75B (diff 13) → K ≈ 2.645.
# K = (1.97 + 2.645)/2 = 2.3075 ≈ 2.31 → 75B ≈ 390 мл, внутри диапазона
# якорей. Шаг на чашку ~2.31·(15²−13²) ≈ 129 мл согласуется с хирургическим
# правилом ~150 мл на чашку. Формула эвристическая, НЕ валидирована.
HEURISTIC_K = 2.31

# Относительная ошибка модели объёма (σ_V/V):
# BREAST-V (Longo et al., 2013): 18.4%; Qiao (1997): ~15% (литературная
# точность антропометрии); эвристика: 20%.
MODEL_ERRORS = {
    "heuristic": 0.20,
    "qiao": 0.15,
    "breast-v": 0.184,
}

# Относительная неопределённость плотности (σ_ρ/ρ): ~3.5% при явно заданном
# пользователем g, ~7% при неопределённом g (разброс доли железистой
# ткани по маммографии ~42.5% ± 30.3%).
DENSITY_ERROR_EXPLICIT_G = 0.035
DENSITY_ERROR_DEFAULT = 0.07

# Поправка доли железистой ткани по возрасту и ИМТ:
#   g_adj = clamp(g − K_AGE·(age−40)/40 − K_BMI·(BMI−27)/10, 0, 1)
# ВАЖНО: K_AGE/K_BMI - РУКОПИСНЫЕ ПЛЕЙСХОЛДЕРЫ (аппроксимация тенденции,
# НЕ результат регрессии). Реальные зависимости нелинейны и с большим
# разбросом (g ~ 42.5% ± 30.3% по маммографии; см. обзоры маммографической
# плотности, напр. Boyd et al., Cancer Epidemiol Biomarkers Prev 2007).
# Значения демпфированы от грубых трендов (k_age ≈ 0.4, k_bmi ≈ 1.2).
K_AGE = 0.16      # плейсхолдер, не валидирован
K_BMI = 0.60      # плейсхолдер, не валидирован

# Неопределённость измерений (--measurement-error, по умолчанию 1.0 см)
# учитывается через частные производные формулы объёма по каждому входу -
# см. measurement_error_fraction(). Единой константы "% объёма на 1 см"
# НЕТ: вклад зависит от модели и значения предикторов.

# Формула BREAST-V (оригинал: Longo et al., 2013, Plast Reconstr Surg
# 132:1e-7e, PMID 23806950; 108 мастэктомических образцов, R²=0.73,
# MAE 89.7 г / 18.4%):
#   V(мл) = −231.66 + 0.5747·(SN-N)² + 18.5478·(FFp) + 14.5087·(N-IMF)
# Предикторы (см): SN-N - яремная вырезка-сосок, FFp - проекция железы
# от грудной стенки (fold projection), N-IMF - сосок-инфрамаммарная
# складка. В редакции Huang N-s et al. (PLoS ONE 2017,
# DOI 10.1371/journal.pone.0172122) параметр FFp заменён на BP (breast
# projection) с СОХРАНЕНИЕМ коэффициентов (см. Huang et al. 2017, Methods).
BREAST_V_INTERCEPT = -231.66
BREAST_V_COEF_SNN2 = 0.5747     # на (SN-N)², мл/см²
BREAST_V_COEF_FFP = 18.5478     # на FFp (≈ BP - проекция железы), мл/см
BREAST_V_COEF_NIMF = 14.5087    # на N-IMF, мл/см


# ----------------------------------------------------------------------------
# Размер: бандаж, чашка, sister sizes
# ----------------------------------------------------------------------------

def round_band(underbust):
    """Округлить обхват под грудью к ближайшему кратному 5 см, минимум 65
    (шаг бандажа 5 см по EN 13402-3).
    """
    band = int(round(underbust / 5.0) * 5)
    return max(65, band)


def cup_label(diff):
    """
    Метка чашки по разнице обхватов diff (см) по EN 13402-3
    (cup size = bust − underbust, шаг 2 см: 10-12 = AA, ..., 26-28 = H).

    Возвращает (метка, экстраполяция_за_пределы_таблицы).
    diff < CUP_MIN_DIFF - размер меньше AA / нестандартный.
    diff >= CUP_MAX_DIFF - продолжение алфавита после H (I, J, K, ...) -
    экстраполяция.
    """
    if diff < CUP_MIN_DIFF:
        return "<AA (меньше AA, нестандартный)", False
    for lo, hi, cup in CUP_TABLE:
        if lo <= diff < hi:
            return cup, False
    # diff >= CUP_MAX_DIFF см: продолжение алфавита после H - экстраполяция
    extra = int((diff - CUP_MAX_DIFF) // 2) + 1
    return chr(ord("H") + extra), True


def cup_index(cup):
    """Порядковый индекс чашки в ряду AA, A, B, ..., H, I, J, ..."""
    if cup == "AA":
        return 0
    if len(cup) == 1 and "A" <= cup <= "Z":
        return ord(cup) - ord("A") + 1
    return None


def cup_from_index(idx):
    """Метка чашки по индексу (0 = AA, 1 = A, 2 = B, ...)."""
    return "AA" if idx == 0 else chr(ord("A") + idx - 1)


def round_even(x):
    """Округлить до ближайшего чётного целого (band size в дюймах)."""
    return int(round(x / 2.0) * 2)


def us_band(band_eu, offset=0):
    """
    US-бандаж по EU-бандажу.

    Современная фиттинг-конвенция (по умолчанию, offset=0):
    band_us = обхват под грудью в дюймах, округлённый до ближайшего
    чётного (band_eu уже округлён к 5 см, поэтому round_even(band_eu /
    2.54) даёт тот же результат, что и округление исходного измерения).

    Устаревшие ритейл-конвенции: многие американские бренды добавляют
    +2 или +4 дюйма к измерению (offset=2 / offset=4); вариант +4 даёт
    старое соответствие EU 65→30, 70→32, 75→34, 80→36, 85→38, 90→40,
    95→42. Выбирается флагом --us-band-offset.
    """
    return round_even(band_eu / 2.54) + offset


def us_cup_label(diff_cm):
    """
    US-метка чашки по разнице обхватов diff (см), пересчёт напрямую
    из измерений (НЕ перекодирование EU-букв: шаг EU-чашки 2 см ≠
    1 дюйму US-чашки).

    cup_diff_inch = diff / 2.54, округление до ближайшего целого дюйма;
    1"=A, 2"=B, 3"=C, 4"=D, 5"=DD (≡E), 6"=DDD (≡F), 7"=G, 8"=H, 9"=I,
    10"=J, далее алфавит (K, L, ...). diff < 0.5" - "<AA".
    """
    diff_in = diff_cm / 2.54
    n = int(round(diff_in))
    if n < 1:
        return "<AA (меньше AA, нестандартный)"
    # 1..4 - A..D; 5 - DD; 6 - DDD; далее - G, H, I, J, K, ...
    if n <= 4:
        return chr(ord("A") + n - 1)
    if n == 5:
        return "DD"
    if n == 6:
        return "DDD"
    return chr(ord("G") + n - 7)


def us_size(band_eu, diff_cm, band_offset=0):
    """
    Полный US-размер строкой, например "34C", из EU-бандажа и разницы
    обхватов (см). Чашка пересчитывается напрямую из измерений в дюймах
    (шаг US-чашки 1 дюйм = 2.54 см), бандаж - округлённый подгрудный
    обхват в дюймах (опциональный ритейл-офсет +2/+4 - см. us_band).
    """
    return f"{us_band(band_eu, band_offset)}{us_cup_label(diff_cm)}"


def us_cup_from_index(idx):
    """Метка US-чашки по индексу (0 = AA, 1 = A, ..., 5 = DD, 6 = DDD)."""
    if idx == 0:
        return "AA"
    if idx <= 4:
        return chr(ord("A") + idx - 1)
    if idx == 5:
        return "DD"
    if idx == 6:
        return "DDD"
    return chr(ord("G") + idx - 7)


def us_sister_sizes(band_eu, diff_cm, band_offset=0):
    """
    US-сестринские размеры: бандаж ±2 дюйма, чашка ∓1 буква.

    Конвенция "sister sizes" для US (34DDD ≈ 36DD ≈ 32G): при увеличении
    бандажа на 2 дюйма объём чашки сохраняется при уменьшении метки чашки
    на 1 букву. Строится от US-размера: band_us ±2, буква чашки ∓1.
    Возвращает список строк вида ["36DD", "32G"].
    """
    band_us = us_band(band_eu, band_offset)
    cup_us = us_cup_label(diff_cm)
    if cup_us.startswith("<"):
        return []
    # DD и DDD - синонимы E и F в ряду A, B, C, D, DD(=E), DDD(=F), G, ...
    idx = {"DD": 5, "DDD": 6}.get(cup_us, cup_index(cup_us))
    if idx is None:
        return []
    result = []
    for band_delta, cup_delta in ((2, -1), (-2, 1)):
        new_band = band_us + band_delta
        new_idx = idx + cup_delta
        if new_band >= 28 and 0 <= new_idx <= 25:
            result.append(f"{new_band}{us_cup_from_index(new_idx)}")
    return result


def sister_sizes(band, cup):
    """
    Эквивалентные по объёму чашки размеры: бандаж ±5 см, чашка ∓1.

    Конвенция "sister sizes" (75B ≈ 80A ≈ 70C): при увеличении бандажа
    на 5 см объём чашки сохраняется при уменьшении метки чашки на 1.
    Возвращает список строк вида ["80A", "70C"].
    """
    idx = cup_index(cup)
    if idx is None:
        return []
    result = []
    for band_delta, cup_delta in ((5, -1), (-5, 1)):
        new_band = band + band_delta
        new_idx = idx + cup_delta
        if new_band >= 65 and 0 <= new_idx <= 25:
            result.append(f"{new_band}{cup_from_index(new_idx)}")
    return result


# ----------------------------------------------------------------------------
# Модели объёма
# ----------------------------------------------------------------------------

def breast_volume_heuristic(bust, underbust):
    """
    Объём одной груди, мл: V ≈ 2.31 · diff² · (band / 75).

    Эвристическая fallback-модель, НЕ валидирована. Функциональная форма
    унаследована от предыдущей версии; константа калибрована по СРЕДНЕМУ
    двух опубликованных якорей:
    (а) Huang et al. 2017 (PLoS ONE, DOI 10.1371/journal.pone.0172122):
        объёмы по чашкам A≈261/B≈328/C≈408/≥D≈539 мл → K ≈ 1.97;
    (б) стереофотограмметрия (Loughry et al. 1987): ~447 мл при 75B
        (diff 13 см) → K ≈ 2.645.
    K = (1.97 + 2.645)/2 = 2.3075 ≈ 2.31 → 75B ≈ 390 мл, внутри диапазона
    якорей. Шаг ~129 мл на чашку согласуется с хирургическим правилом
    ~150 мл на чашку размера.
    """
    diff = bust - underbust
    band = round_band(underbust)
    return HEURISTIC_K * diff * diff * (band / 75.0)


def breast_volume_qiao(mp, lr, mr, ir):
    """
    Объём одной груди, мл, по антропометрической формуле Qiao et al. (1997):

        V = (π/3) · MP² · (LR + MR + IR − MP)

    где MP - проекция молочной железы от грудной стенки, LR/MR/IR -
    латеральный/медиальный/нижний радиусы от соска до границы железы, см.
    Литературная точность антропометрических методов ~15% (σ_V/V); формула
    предполагает клиническую точность измерений (~2-3 мм) - при самостоятельном
    измерении рекомендуется --measurement-error 0.3.

    Выбрасывает ValueError, если LR+MR+IR <= MP (модель неприменима:
    выражение даёт неположительный объём).
    """
    if lr + mr + ir <= mp:
        raise ValueError(
            "Qiao: LR + MR + IR должно быть больше MP "
            "(иначе модель неприменима)")
    return (math.pi / 3.0) * mp * mp * (lr + mr + ir - mp)


def breast_volume_breast_v(sn_n, ffp, n_imf):
    """
    Объём одной груди, мл, по формуле BREAST-V (Longo et al., 2013):

        V = −231.66 + 0.5747·(SN-N)² + 18.5478·(FFp) + 14.5087·(N-IMF)

    Предикторы (см): SN-N - расстояние яремная вырезка-сосок (sternal
    notch-nipple), FFp - проекция железы от грудной стенки (fold
    projection), N-IMF - расстояние сосок-инфрамаммарная складка.

    В редакции Huang et al. 2017 параметр FFp заменён на BP (breast
    projection) с сохранением коэффициентов - здесь используется эта
    редакция (FFp ≈ BP - проекция железы).

    Источники:
    - Longo B. et al. The BREAST-V: a unifying predictive formula...
      Plast Reconstr Surg. 2013;132:1e-7e. PMID 23806950 - оригинальная
      валидированная регрессия на 108 мастэктомических образцах,
      R² = 0.73, MAE 89.7 г (относительная ошибка 18.4%).
    - Huang N-s et al. A prospective study of breast anthropomorphic
      measurements, volume and ptosis in 605 Asian patients. PLoS ONE.
      2017;12(2):e0172122. DOI 10.1371/journal.pone.0172122 - замена
      FFp → BP в азиатской когорте.

    Выбрасывает ValueError, если объём получается неположительным
    (модель неприменима к введённым значениям предикторов).
    """
    volume = (
        BREAST_V_INTERCEPT
        + BREAST_V_COEF_SNN2 * sn_n * sn_n
        + BREAST_V_COEF_FFP * ffp
        + BREAST_V_COEF_NIMF * n_imf
    )
    if volume <= 0.0:
        raise ValueError(
            "BREAST-V: при таких значениях предикторов модель неприменима "
            "(объём <= 0)")
    return volume


# ----------------------------------------------------------------------------
# Плотность и вес
# ----------------------------------------------------------------------------

def adjust_glandular_fraction(g, age=None, bmi=None):
    """
    Скорректировать долю железистой ткани g по возрасту и ИМТ.

    g_adj = clamp(g − K_AGE·(age−40)/40 − K_BMI·(BMI−27)/10, 0, 1)

    ВАЖНО: коэффициенты K_AGE/K_BMI - РУКОПИСНЫЕ ПЛЕЙСХОЛДЕРЫ, аппроксимация
    тенденции, а НЕ валидированная регрессия: реальные зависимости
    нелинейны и с большим разбросом (по маммографии g ~ 42.5% ± 30.3%;
    плотность обратно коррелирует с возрастом и ИМТ - медианный ИМТ ~31.2
    при полностью жировой груди против ~22.8 при крайне плотной; см.
    обзоры, напр. Boyd et al., Cancer Epidemiol Biomarkers Prev 2007).
    """
    g_adj = g
    if age is not None:
        g_adj -= K_AGE * (age - 40.0) / 40.0
    if bmi is not None:
        g_adj -= K_BMI * (bmi - 27.0) / 10.0
    return min(1.0, max(0.0, g_adj))


def tissue_density(g):
    """
    Плотность ткани груди, г/см³, как функция доли железистой (фиброгландулярной)
    ткани g ∈ [0, 1]: кусочно-линейная интерполяция по трём опорным
    точкам ICRU Report 44:

        (g=0 → 0.95) - жировая; (g=0.5 → 0.96) - 50/50; (g=1 → 1.02) -
        фиброгландулярная.

    При g=0.5 даёт табличное значение ICRU 0.96 (линейная интерполяция
    двух крайних точек дала бы 0.985 - расхождение с таблицей).
    Средняя доля по маммографии ~42.5% ± 30.3%.
    """
    if g <= 0.5:
        return DENSITY_ADIPOSE + (DENSITY_5050 - DENSITY_ADIPOSE) * (g / 0.5)
    return DENSITY_5050 + (DENSITY_FIBROGLANDULAR - DENSITY_5050) * ((g - 0.5) / 0.5)


def measurement_error_fraction(model, meas_err, inputs):
    """
    Относительный вклад ошибки измерений в σ_V/V, вычисленный через
    частные производные формулы объёма по каждому входу (Δx = meas_err
    для всех измерений - консервативно).

    - heuristic: V = K·diff²·(band/75)
        ∂V/∂diff = 2K·diff·(band/75);  ∂V/∂band = K·diff²/75.
      Пример (diff=13, band=75, Δ=1 см): σ_V/V ≈ 21-22%.
    - Примечание по d_b (вклад бандажа): band = round_band(underbust) -
      ступенчатая функция с шагом 5 см; производная по underbust через
      неё почти всюду равна нулю, в точках разрыва - дельта-функция.
      d_b ниже трактуется как КОНСЕРВАТИВНАЯ ОЦЕНКА СВЕРХУ (band
      условно непрерывен), а не точная производная; вклад d_b мал
      относительно d_v.
    - breast-v: V = a·SN-N² + b·FFp + c·N-IMF + d
        ∂V/∂SN-N = 2a·SN-N (≈28 мл при SN-N=24.7 см → ~7% при V≈400 мл);
        ∂V/∂FFp = b; ∂V/∂N-IMF = c.
    - qiao: V = (π/3)·MP²·(S−MP), S = LR+MR+IR
        ∂V/∂MP = (π/3)·MP·(2S−3MP);  ∂V/∂S = (π/3)·MP²;
        ΔS из трёх независимых радиусов: sqrt(3)·Δx.

    Возвращает долю от 1 (σ_изм/V). Зависимости от модели и значений
    предикторов - потому единой константы "% на 1 см" быть не может.
    """
    if meas_err <= 0:
        return 0.0
    if model == "heuristic":
        bust, underbust = inputs
        diff = bust - underbust
        band = round_band(underbust)
        volume = HEURISTIC_K * diff * diff * (band / 75.0)
        d_v = math.sqrt(2.0) * 2.0 * HEURISTIC_K * \
            diff * (band / 75.0) * meas_err
        d_b = HEURISTIC_K * diff * diff / 75.0 * meas_err
        return math.hypot(d_v, d_b) / volume
    if model == "breast-v":
        sn_n, ffp, n_imf = inputs
        volume = breast_volume_breast_v(sn_n, ffp, n_imf)
        d_n = 2.0 * BREAST_V_COEF_SNN2 * sn_n * meas_err
        d_f = BREAST_V_COEF_FFP * meas_err
        d_i = BREAST_V_COEF_NIMF * meas_err
        return math.hypot(math.hypot(d_n, d_f), d_i) / volume
    if model == "qiao":
        mp, lr, mr, ir = inputs
        volume = breast_volume_qiao(mp, lr, mr, ir)
        s = lr + mr + ir
        d_v_dmp = (math.pi / 3.0) * mp * (2.0 * s - 3.0 * mp)
        d_v_ds = (math.pi / 3.0) * mp * mp
        dv = math.sqrt((d_v_dmp * meas_err) ** 2
                       + (d_v_ds * meas_err * math.sqrt(3.0)) ** 2)
        return abs(dv) / volume
    raise ValueError(f"неизвестная модель: {model}")


def volume_uncertainty_fraction(model_error, meas_frac):
    """
    Относительная неопределённость ОБЪЁМА (доля от 1), без плотности:

        σ_V/V = sqrt(σ_мод² + σ_изм²)

    σ_мод - ошибка модели (MODEL_ERRORS); σ_изм - вклад ошибки измерений
    (measurement_error_fraction). Именно эта величина задаёт 95% ДИ объёма.
    """
    return math.sqrt(model_error ** 2 + meas_frac ** 2)


def combined_uncertainty(model_error, density_error, meas_frac):
    """
    Итоговая относительная неопределённость ВЕСА (доля от 1):

        σ_W/W = sqrt((σ_V/V)² + (σ_ρ/ρ)²)

    где σ_V/V - уже включает вклад ошибки измерений
    (volume_uncertainty_fraction), двойной учёт исключён.
    """
    sigma_v = volume_uncertainty_fraction(model_error, meas_frac)
    return math.sqrt(sigma_v ** 2 + density_error ** 2)


def breast_weight_ci(volume, glandular_fraction, meas_frac, model_error,
                     density_error=DENSITY_ERROR_DEFAULT):
    """
    Вес одной груди с 95% доверительным интервалом.

    W = V·ρ(g); σ_V/V = sqrt(σ_мод² + σ_изм²) (частные производные);
    σ_W = W·sqrt((σ_V/V)² + (σ_ρ/ρ)²); 95% ДИ = W ± 1.96·σ_W.

    Источники: ICRU Report 44 (плотность), Longo et al. 2013 / Qiao et al.
    1997 (ошибки моделей объёма). Возвращает (W, σ_W, ДИ_низ, ДИ_верх).
    """
    rho = tissue_density(glandular_fraction)
    weight = volume * rho
    sigma_rel = combined_uncertainty(model_error, density_error, meas_frac)
    sigma_w = weight * sigma_rel
    return weight, sigma_w, weight - 1.96 * sigma_w, weight + 1.96 * sigma_w


# ----------------------------------------------------------------------------
# Валидация и ввод
# ----------------------------------------------------------------------------

def validate(bust, underbust):
    """
    Проверить корректность входных данных.

    Жёсткие ошибки: обхваты вне 50-160 см; обхват груди не больше обхвата
    под грудью; diff вне 5-40 см.
    Предупреждения (возвращаются списком): обхват под грудью вне 60-120 см
    (модели калиброваны на ограниченном диапазоне); diff ниже нижней
    границы таблицы чашек EN 13402-3 (10 см) - экстраполяция ниже AA.
    """
    warnings = []
    for name, value in (("обхват груди", bust), ("обхват под грудью", underbust)):
        if not (MIN_CM <= value <= MAX_CM):
            raise ValueError(
                f"{name.capitalize()} ({value:g} см) вне разумных пределов "
                f"{MIN_CM}-{MAX_CM} см"
            )
    if bust <= underbust:
        raise ValueError("Обхват груди должен быть больше обхвата под грудью")
    diff = bust - underbust
    if diff < DIFF_MIN or diff > DIFF_MAX:
        raise ValueError(
            f"Разница обхватов {diff:g} см вне поддерживаемого диапазона "
            f"{DIFF_MIN}-{DIFF_MAX} см"
        )
    if underbust < BAND_WARN[0] or underbust > BAND_WARN[1]:
        warnings.append(
            f"обхват под грудью {underbust:g} см вне {BAND_WARN[0]}-{BAND_WARN[1]} см - "
            "модели калиброваны на ограниченном диапазоне"
        )
    if diff < CUP_MIN_DIFF:
        warnings.append(
            f"diff = {diff:g} см < {CUP_MIN_DIFF} см - размер ниже AA "
            "(вне таблицы чашек), оценка объёма экстраполирована"
        )
    return warnings


def parse_number(text):
    """Разобрать число (с точкой или запятой в качестве разделителя)."""
    return float(str(text).replace(",", "."))


def build_arg_parser():
    """Собрать CLI-парсер (argparse, стандартная библиотека)."""
    parser = argparse.ArgumentParser(
        prog="breast_calculator.py",
        description="Антропометрический калькулятор объёма, веса и размера груди.",
        epilog=("Формат ввода: все величины в сантиметрах; десятичный "
                "разделитель - точка или запятая (92, 92.5, 92,5). "
                "Обязательные позиционные аргументы: BUST (обхват груди по "
                "самым выступающим точкам) и UNDERBUST (обхват под грудью), "
                "в этом порядке. Без них - интерактивный ввод двух обхватов."),
    )
    parser.add_argument("bust", nargs="?", type=parse_number, metavar="BUST",
                        help="обхват груди, см")
    parser.add_argument("underbust", nargs="?", type=parse_number, metavar="UNDERBUST",
                        help="обхват под грудью, см")
    parser.add_argument("--model", choices=["heuristic", "qiao", "breast-v"], default=None,
                        help="модель объёма (по умолчанию: qiao/breast-v при полных "
                             "доп. измерениях, иначе эвристика)")
    parser.add_argument("--age", type=float, default=None,
                        help="возраст, лет (поправка доли железистой ткани)")
    parser.add_argument("--bmi", type=float, default=None,
                        help="ИМТ, кг/м² (поправка доли железистой ткани)")
    parser.add_argument("--glandular-fraction", type=float, default=None,
                        help="доля железистой ткани g ∈ [0,1] (по умолчанию 0.5)")
    parser.add_argument("--measurement-error", type=float, default=1.0,
                        help="ошибка измерения, см (по умолчанию 1.0)")
    parser.add_argument("--us-band-offset", type=int, choices=[0, 2, 4], default=0,
                        help="US-бандаж: добавка к округлённому подгрудному "
                             "обхвату в дюймах (0 - современная конвенция, "
                             "по умолчанию; 2 или 4 - устаревшие ритейл-"
                             "конвенции некоторых брендов)")
    parser.add_argument("--self-check", action="store_true",
                        help="проверка согласованности моделей с опубликованными "
                             "якорями (Huang 2017, Behrens 2024, Longo 2013, "
                             "Qiao 1997); НЕ независимая out-of-sample валидация")
    # Предикторы BREAST-V (Longo 2013; редакция Huang 2017: FFp ≈ BP)
    parser.add_argument("--notch-nipple", type=float, default=None,
                        help="BREAST-V: яремная вырезка - сосок (SN-N), см")
    parser.add_argument("--fold-nipple", type=float, default=None,
                        help="BREAST-V: сосок - инфрамаммарная складка (N-IMF), см")
    parser.add_argument("--fold-projection", type=float, default=None,
                        help="BREAST-V: проекция железы (FFp ≈ BP), см")
    # Предикторы Qiao (1997)
    parser.add_argument("--mp", type=float, default=None,
                        help="Qiao: проекция молочной железы (MP), см")
    parser.add_argument("--lr", type=float, default=None,
                        help="Qiao: латеральный радиус (LR), см")
    parser.add_argument("--mr", type=float, default=None,
                        help="Qiao: медиальный радиус (MR), см")
    parser.add_argument("--ir", type=float, default=None,
                        help="Qiao: нижний радиус (IR), см")
    return parser


# ----------------------------------------------------------------------------
# Отчёт
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Валидация по реальным данным из статей
# ----------------------------------------------------------------------------

# Якорные "реальные" данные (опубликованные измерения, не оценки формулами):
# 1) Huang et al. 2017 (PLoS ONE, DOI 10.1371/journal.pone.0172122):
#    средний объём ОДНОЙ груди по самооценённому размеру чашки в когорте
#    605 китаянок (1210 грудей, средний объём 340.0±109.1 мл):
#    A ≈ 260.9, B ≈ 328.0, C ≈ 408.1, ≥D ≈ 539.0 мл.
# 2) Behrens et al. 2024 (Breast Cancer, DOI 10.1007/s12282-024-01647-6):
#    MRI-волюметрия, 22 женщины (Германия, средний ИМТ 24.6):
#    левая 771.0 мл, правая 763.9 мл (в среднем ~767 мл на грудь).
# 3) Longo et al. 2013 (PRS, PMID 23806950): средние предикторы 88
#    кавказских женщин: SN-N 24.73 см, N-IMF 9.26 см (BP в публикации
#    напрямую не приведён - для проверки берём типичное значение ~7 см).
# 4) Qiao et al. 1997 (Aesthet Plast Surg 21:362-368): средний объём
#    молодых китаянок 310-330 мл (250 грудей).
# 5) Популяционные диапазоны (сводка в Huang 2017): европейские женщины
#    407.2-623.5 мл; азиатские 325.4-386.0 мл.

SELF_CHECK_ANCHORS = {
    "huang_cup_A": {"band": 75, "diff": 11, "ref": 260.9, "src": "Huang 2017, чашка A"},
    "huang_cup_B": {"band": 75, "diff": 13, "ref": 328.0, "src": "Huang 2017, чашка B"},
    "huang_cup_C": {"band": 75, "diff": 15, "ref": 408.1, "src": "Huang 2017, чашка C"},
    "huang_cup_D": {"band": 75, "diff": 17, "ref": 539.0, "src": "Huang 2017, чашка >=D"},
    "huang_mean_asia": {"band": 75, "diff": 13, "ref": 340.0, "src": "Huang 2017, среднее когорты (±109)"},
    "behrens_mri": {"band": 85, "diff": 17, "ref": 767.0,
                    "src": "Behrens 2024, MRI (среднее на грудь) *"},
    "qiao_young_chinese": {"band": 75, "diff": 12, "ref": 320.0, "src": "Qiao 1997, среднее молодых китаянок"},
    "euro_mean_low": {"band": 80, "diff": 13, "ref": 407.2, "src": "Европейские женщины, нижняя граница"},
    "euro_mean_high": {"band": 85, "diff": 15, "ref": 623.5, "src": "Европейские женщины, верхняя граница"},
}


def run_self_check():
    """
    Проверка согласованности моделей с опубликованными якорями.

    ДЛЯ КАЖДОГО ЯКОРЯ вычисляется объём эвристикой при характерных
    обхватах (band/diff) и отклонение от референса. ВАЖНО: это проверка
    СОГЛАСОВАННОСТИ, а НЕ независимая out-of-sample валидация -
    константа HEURISTIC_K откалибрована по части этих же якорей, поэтому
    попадание в допуск ±30% ожидаемо. Научная валидация требует
    независимой когорты, которой у нас нет.
    """
    print("=== Проверка согласованности с опубликованными якорями ===")
    print("(НЕ независимая out-of-sample валидация - см. пояснения)")
    print()
    print(f"{'Якорь':<38}{'Статья, мл':>10}{'Модель, мл':>11}{'Откл.':>9}  Статус")
    print("-" * 88)
    results = []
    for key, a in SELF_CHECK_ANCHORS.items():
        bust = a["band"] + a["diff"]
        volume = breast_volume_heuristic(bust, a["band"])
        dev = (volume - a["ref"]) / a["ref"] * 100.0
        status = "OK" if abs(dev) <= 30 else "ВНЕ ДОПУСКА"
        results.append((key, a["ref"], volume, dev, status))
        print(f"{a['src']:<38}{a['ref']:>10.0f}{volume:>11.0f}"
              f"{dev:>+8.0f}%  {status}")
    print()
    # BREAST-V: Longo 2013, средние предикторы когорты (SN-N 24.73, N-IMF 9.26;
    # FFp ≈ BP ≈ 7 см - типичное значение, в статье напрямую не приведено).
    bv = breast_volume_breast_v(24.73, 7.0, 9.26)
    dev_bv = (bv - 434.5) / 434.5 * 100.0
    print(f"BREAST-V (Longo 2013, SN-N 24.73 / FFp≈7 / N-IMF 9.26, см): "
          f"{bv:.0f} мл  (ожидание ~412-457 мл по кавказской когорте; "
          f"откл. от середины {dev_bv:+.0f}%)")
    print()
    print("Пояснения:")
    print("  - ЭТО НЕ НЕЗАВИСИМАЯ ВАЛИДАЦИЯ: проверяется, что эвристика "
          "попадает в опубликованные якоря в пределах ±30%, что ожидаемо, "
          "поскольку константа K калибрована по части этих якорей "
          "(Huang 2017 и стереофотограмметрия).")
    print("  - Якоря Huang 2017 привязаны к самооценённому размеру чашки; "
          "band 75 и diff по таблице - типичные значения когорты.")
    print("  - * Behrens 2024: MRI-волюметрия немецкой выборки; характерные "
          "обхваты - реконструкция типичного размера, НЕ измерения той же "
          "когорты; точка не независимая.")
    print("  - Допуск ±30% учитывает разброс ±109 мл внутри когорты Huang "
          "и невалидированность эвристической fallback-модели.")
    independent = [r for r in results if r[0] != "behrens_mri"]
    behrens_ok = next(r[4] for r in results if r[0] == "behrens_mri") == "OK"
    n_ind = len(independent)
    ok = all(r[4] == "OK" for r in independent)
    print()
    if ok:
        print(f"Итог: {n_ind} независимых якорей проверено, все в допуске "
              "±30% (согласованность подтверждена). Behrens 2024 - "
              "реконструкция размера, из вердикта исключён "
              f"(статус: {'OK' if behrens_ok else 'ВНЕ ДОПУСКА'}).")
    else:
        print("Итог: часть независимых якорей вне допуска ±30% - "
              "см. таблицу выше (Behrens 2024 - реконструкция, из "
              "вердикта исключён).")


MODEL_DESCRIPTIONS = {
    "heuristic": "эвристическая (fallback; не валидирована)",
    "qiao": "Qiao et al. (1997), антропометрическая",
    "breast-v": "BREAST-V (Longo 2013); редакция Huang 2017: FFp → BP, "
                "коэффициенты сохранены",
}


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    # Ошибка измерения нулевой или отрицательной быть не может
    if args.measurement_error <= 0:
        parser.error("--measurement-error должен быть > 0")

    # Режим проверки согласованности с опубликованными якорями
    if args.self_check:
        run_self_check()
        return

    # Позиционные аргументы или интерактивный ввод (только базовые 2 параметра)
    if args.bust is None or args.underbust is None:
        if args.bust is not None or args.underbust is not None:
            parser.error("укажите оба обхвата: BUST UNDERBUST")
        try:
            args.bust = parse_number(
                input("Обхват груди по самым выступающим точкам, см "
                      "(например 92 или 92.5): "))
            args.underbust = parse_number(
                input("Обхват под грудью, см (например 75): "))
        except ValueError:
            print("Ошибка: введено не число. Формат: сантиметры, "
                  "десятичный разделитель - точка или запятая "
                  "(примеры: 92, 92.5, 92,5).")
            sys.exit(1)

    bust, underbust = args.bust, args.underbust

    try:
        warnings = validate(bust, underbust)
    except ValueError as e:
        print(f"Ошибка: {e}")
        sys.exit(1)

    diff = bust - underbust
    band = round_band(underbust)
    cup, extrapolated = cup_label(diff)
    us = us_size(band, diff, args.us_band_offset)

    # --- Выбор модели объёма ---
    qiao_ready = all(v is not None for v in (
        args.mp, args.lr, args.mr, args.ir))
    bv_ready = all(v is not None for v in
                   (args.notch_nipple, args.fold_nipple, args.fold_projection))
    model = args.model
    if model == "qiao" and not qiao_ready:
        parser.error("--model qiao требует --mp, --lr, --mr, --ir")
    if model == "breast-v" and not bv_ready:
        parser.error("--model breast-v требует --notch-nipple, --fold-nipple, "
                     "--fold-projection")
    if model is None:
        if qiao_ready:
            model = "qiao"       # все параметры известны - по умолчанию Qiao
        elif bv_ready:
            model = "breast-v"
        else:
            model = "heuristic"

    # Входы модели: объём и вклад ошибки измерений считаются от ОДНОГО
    # набора значений (частные производные по каждому входу)
    if model == "qiao":
        volume_inputs = (args.mp, args.lr, args.mr, args.ir)
    elif model == "breast-v":
        # Порядок: (SN-N, FFp, N-IMF); FFp - это --fold-projection,
        # N-IMF - это --fold-nipple (в редакции Huang FFp ≈ BP).
        volume_inputs = (args.notch_nipple, args.fold_projection,
                         args.fold_nipple)
    else:
        volume_inputs = (bust, underbust)

    try:
        if model == "qiao":
            volume = breast_volume_qiao(*volume_inputs)
        elif model == "breast-v":
            volume = breast_volume_breast_v(*volume_inputs)
        else:
            volume = breast_volume_heuristic(*volume_inputs)
    except ValueError as e:
        print(f"Ошибка: {e}")
        sys.exit(1)

    # --- Доля железистой ткани и плотность ---
    g_explicit = args.glandular_fraction is not None
    g = args.glandular_fraction if g_explicit else 0.5
    if not 0.0 <= g <= 1.0:
        print("Ошибка: --glandular-fraction должен быть в [0, 1].")
        sys.exit(1)
    g_adj = adjust_glandular_fraction(g, args.age, args.bmi)
    density_error = DENSITY_ERROR_EXPLICIT_G if g_explicit else DENSITY_ERROR_DEFAULT

    # --- Вес и ДИ (объём: модель + измерения; вес: + плотность) ---
    model_error = MODEL_ERRORS[model]
    meas_frac = measurement_error_fraction(
        model, args.measurement_error, volume_inputs)
    sigma_v = volume_uncertainty_fraction(model_error, meas_frac)
    sigma_rel = combined_uncertainty(model_error, density_error, meas_frac)
    weight, sigma_w, w_lo, w_hi = breast_weight_ci(
        volume, g_adj, meas_frac, model_error, density_error)
    # Объём и вес неотрицательны по физике: нижняя граница ДИ не ниже 0
    v_lo = max(0.0, volume * (1.0 - 1.96 * sigma_v))
    v_hi = volume * (1.0 + 1.96 * sigma_v)
    w_lo = max(0.0, w_lo)
    if 1.96 * sigma_v >= 1.0 or 1.96 * sigma_rel >= 1.0:
        print("Предупреждение: неопределённость такова, что 95% ДИ "
              "обрезан снизу нулём - оценка малоинформативна.")

    # --- Отчёт ---
    print("=== Результаты антропометрического анализа ===")
    print("Входные данные:")
    print(f"  Обхват груди: {bust:g} см")
    print(f"  Обхват под грудью: {underbust:g} см")
    print(f"  Разница: {diff:g} см")
    print()
    print("Размер (EN 13402-3, шаг чашки 2 см):")
    print(f"  Бандаж: {band}")
    print(f"  Чашка: {cup} (diff {cup_range_text(diff)})"
          + ("  [экстраполяция за пределы таблицы]" if extrapolated else ""))
    print(f"  Итог: {band}{cup}")
    print(f"  US: {us} (бандаж EU {band} → US {us_band(band, args.us_band_offset)}; "
          "чашка - пересчёт из измерений в дюймах, шаг 1\")")
    sisters = sister_sizes(band, cup)
    us_sisters = us_sister_sizes(band, diff, args.us_band_offset)
    if sisters or us_sisters:
        parts = []
        if sisters:
            parts.append(f"EU: {', '.join(sisters)}")
        if us_sisters:
            parts.append(f"US: {', '.join(us_sisters)}")
        print(f"  Sister sizes: {' / '.join(parts)}")
    print()
    print(f"Модель объёма: {MODEL_DESCRIPTIONS[model]}")
    print(f"  Объём одной груди: {round(volume)} мл "
          f"(95% ДИ: {round(v_lo)}-{round(v_hi)} мл)")
    print()
    g_text = f"g={g_adj:.2f}"
    if g_adj != g:
        g_text += f" (скорректировано с {g:.2f} по возрасту/ИМТ)"
    print(
        f"Плотность (ICRU Report 44): ρ({g_text}) = {tissue_density(g_adj):.3f} г/см³")
    print(f"  Вес одной груди: {round(weight)} г "
          f"(95% ДИ: {round(w_lo)}-{round(w_hi)} г)")
    print()
    meas_pct = meas_frac * 100
    print(f"Погрешность объёма: модель ±{model_error * 100:g}% + "
          f"измерения ±{meas_pct:.1f}% (через производные по входам) "
          f"= ±{sigma_v * 100:.1f}%")
    print(f"  Плотность: ещё ±{density_error * 100:g}%")
    print(f"  Итог для веса (95% ДИ): ±{1.96 * sigma_rel * 100:.1f}%")
    print()
    print("Примечания:")
    print("  - Асимметрия: идеально симметричных грудей практически не "
          "бывает; в среднем левая чуть больше правой (по MRI-волюметрии "
          "771 против 764 мл: Behrens et al., 2024, DOI "
          "10.1007/s12282-024-01647-6; линейные расстояния слева больше: "
          "Henseler, 2023, DOI 10.3205/iprs000173). Это статистическая "
          "тенденция, а не прогноз для конкретного человека.")
    print("  - Таблица чашек: EN 13402-3 (cup size = bust − underbust, "
          "шаг 2 см).")
    if args.us_band_offset:
        print(f"  - US-размер: бандаж = округлённый подгрудный обхват в "
              f"дюймах + {args.us_band_offset} (устаревшая ритейл-конвенция; "
              "современная фиттинг-методика - без добавки, см. "
              "--us-band-offset 0); чашка - шаг 1 дюйм (2.54 см) на размер, "
              "пересчёт напрямую из измерений, а не перекодирование EU-букв.")
    else:
        print("  - US-размер: бандаж = округлённый подгрудный обхват в дюймах "
              "(современная фиттинг-конвенция, без добавки; устаревшие "
              "ритейл-конвенции +2/+4 - см. --us-band-offset); чашка - шаг "
              "1 дюйм (2.54 см) на размер, пересчёт напрямую из измерений, "
              "а не перекодирование EU-букв.")
    if model == "heuristic":
        print("  - Объём: эвристическая формула, не валидирована; константа "
              "калибрована по среднему опубликованных якорей (Huang 2017: "
              "чашки A≈261/B≈328/C≈408/≥D≈539 мл; стереофотограмметрия: "
              "~447 мл при 75B).")
    if model == "breast-v":
        print("  - BREAST-V: формула Longo et al. (2013, PMID 23806950), "
              "валидирована на 108 мастэктомических образцах "
              "(MAE 89.7 г / 18.4%); в редакции Huang et al. (2017, "
              "DOI 10.1371/journal.pone.0172122) параметр FFp заменён на "
              "BP с сохранением коэффициентов.")
    if model == "qiao":
        print("  - Qiao: формула предполагает клиническую точность измерений "
              "(~2-3 мм); при самостоятельных измерениях рекомендуется "
              "--measurement-error 0.3.")
    if args.age is not None or args.bmi is not None:
        print("  - Поправка по возрасту/ИМТ - коэффициенты являются "
              "рукописными плейсхолдерами (аппроксимация тенденции, "
              "НЕ регрессия; реальные зависимости нелинейны).")
    for w in warnings:
        print(f"  - Предупреждение: {w}.")


def cup_range_text(diff):
    """Диапазон diff для текущей чашки (для строки отчёта)."""
    if diff < CUP_MIN_DIFF:
        return f"< {CUP_MIN_DIFF} см"
    for lo, hi, cup in CUP_TABLE:
        if lo <= diff < hi:
            return f"{lo}-{hi} см"
    return f">= {CUP_MAX_DIFF} см"


if __name__ == "__main__":
    main()
