# app.py - Streamlit Cloud учун мослаштирилган генетик синдромлар хавф бахолаш дастури
# DELFIA Revvity реагентлари асосида
# Даун, Эдвардс, Патау, Тернер, НТД ва ёш хавфлари учун тулик дастур

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from datetime import datetime, date
import json
import base64
import warnings
import os
warnings.filterwarnings('ignore')

# Streamlit Cloud конфигурацияси
__import__('pysqlite3')
import sys
sys.modules['sqlite3'] = sys.modules.pop('pysqlite3')

# ==================== ЎЗГАРМАСЛАР ====================

# Генетик синдромлар учун асосий хавфлар
BASE_RISKS = {
    'downs': 1/800,      # Даун синдроми
    'edwards': 1/3000,   # Эдвардс синдроми
    'patau': 1/5000,     # Патау синдроми
    'turner': 1/2500,    # Тернер синдроми
    'ntd': 1/1000        # Нейротубуляр дефект
}

# Ёш бўйича хавф кўпайтирувчилари
AGE_MULTIPLIERS = {
    20: {'downs': 0.5, 'edwards': 0.3, 'patau': 0.3, 'turner': 0.4},
    25: {'downs': 0.7, 'edwards': 0.5, 'patau': 0.5, 'turner': 0.6},
    30: {'downs': 1.0, 'edwards': 1.0, 'patau': 1.0, 'turner': 1.0},
    35: {'downs': 2.5, 'edwards': 3.0, 'patau': 3.5, 'turner': 2.0},
    40: {'downs': 5.0, 'edwards': 8.0, 'patau': 10.0, 'turner': 4.0},
    45: {'downs': 10.0, 'edwards': 15.0, 'patau': 20.0, 'turner': 8.0}
}

# DELFIA Revvity нормалари (хақиқий референс қийматлари)
DELFIA_FIRST_TRIMESTER = {
    'PAPP_A': {
        'unit': 'U/L',
        'ranges_by_week': {
            10: {'min': 0.4, 'max': 3.0, 'median': 1.0},
            11: {'min': 0.5, 'max': 3.5, 'median': 1.2},
            12: {'min': 0.6, 'max': 4.0, 'median': 1.4},
            13: {'min': 0.7, 'max': 4.5, 'median': 1.6},
            14: {'min': 0.8, 'max': 5.0, 'median': 1.8}
        },
        'MoM_range': {'low': 0.4, 'high': 2.5}
    },
    
    'FREE_BETA_HCG': {
        'unit': 'ng/ml',
        'ranges_by_week': {
            10: {'min': 15.0, 'max': 120.0, 'median': 40.0},
            11: {'min': 20.0, 'max': 150.0, 'median': 60.0},
            12: {'min': 25.0, 'max': 180.0, 'median': 80.0},
            13: {'min': 30.0, 'max': 200.0, 'median': 100.0},
            14: {'min': 35.0, 'max': 220.0, 'median': 120.0}
        },
        'MoM_range': {'low': 0.5, 'high': 2.0}
    },
    
    'NT': {
        'unit': 'мм',
        'ranges_by_week': {
            10: {'min': 0.8, 'max': 2.2, 'median': 1.2},
            11: {'min': 0.8, 'max': 2.5, 'median': 1.3},
            12: {'min': 0.8, 'max': 2.8, 'median': 1.4},
            13: {'min': 0.8, 'max': 3.0, 'median': 1.5},
            14: {'min': 0.8, 'max': 3.0, 'median': 1.5}
        },
        'normal_max': 2.5
    }
}

DELFIA_SECOND_TRIMESTER = {
    'AFP': {
        'unit': 'ng/ml',
        'ranges_by_week': {
            15: {'min': 15.0, 'max': 60.0, 'median': 30.0},
            16: {'min': 17.0, 'max': 65.0, 'median': 35.0},
            17: {'min': 20.0, 'max': 70.0, 'median': 40.0},
            18: {'min': 22.0, 'max': 75.0, 'median': 45.0},
            19: {'min': 25.0, 'max': 80.0, 'median': 50.0},
            20: {'min': 27.0, 'max': 85.0, 'median': 55.0},
            21: {'min': 30.0, 'max': 90.0, 'median': 60.0},
            22: {'min': 32.0, 'max': 95.0, 'median': 65.0}
        },
        'MoM_range': {'low': 0.5, 'high': 2.0}
    },
    
    'TOTAL_HCG': {
        'unit': 'IU/L',
        'ranges_by_week': {
            15: {'min': 10000, 'max': 60000, 'median': 30000},
            16: {'min': 8000, 'max': 55000, 'median': 28000},
            17: {'min': 7000, 'max': 50000, 'median': 25000},
            18: {'min': 6000, 'max': 45000, 'median': 22000},
            19: {'min': 5000, 'max': 40000, 'median': 20000},
            20: {'min': 4000, 'max': 35000, 'median': 18000},
            21: {'min': 3500, 'max': 30000, 'median': 16000},
            22: {'min': 3000, 'max': 25000, 'median': 14000}
        },
        'MoM_range': {'low': 0.5, 'high': 2.0}
    },
    
    'UE3': {
        'unit': 'nmol/L',
        'ranges_by_week': {
            15: {'min': 1.0, 'max': 5.0, 'median': 2.5},
            16: {'min': 1.5, 'max': 6.0, 'median': 3.0},
            17: {'min': 2.0, 'max': 7.0, 'median': 3.5},
            18: {'min': 2.5, 'max': 8.0, 'median': 4.0},
            19: {'min': 3.0, 'max': 9.0, 'median': 4.5},
            20: {'min': 3.5, 'max': 10.0, 'median': 5.0},
            21: {'min': 4.0, 'max': 11.0, 'median': 5.5},
            22: {'min': 4.5, 'max': 12.0, 'median': 6.0}
        },
        'MoM_range': {'low': 0.5, 'high': 2.0}
    }
}

# ==================== ФУНКЦИЯЛАР ====================

def calculate_bmi(weight, height):
    """BMI ҳисоблаш"""
    if height > 0:
        return round(weight / ((height/100) ** 2), 1)
    return 22.0

def get_delfia_norm(parameter, gestational_week, trimester="first"):
    """DELFIA Revvity нормаларини олиш"""
    if trimester == "first":
        norms = DELFIA_FIRST_TRIMESTER
    else:
        norms = DELFIA_SECOND_TRIMESTER
    
    if parameter in norms and gestational_week in norms[parameter]['ranges_by_week']:
        return norms[parameter]['ranges_by_week'][gestational_week]
    else:
        closest_week = min(norms[parameter]['ranges_by_week'].keys(), 
                          key=lambda x: abs(x - gestational_week))
        return norms[parameter]['ranges_by_week'][closest_week]

def calculate_mom_delfia(value, parameter, gestational_week, maternal_weight=None, trimester="first"):
    """DELFIA Revvity учун MoM ҳисоблаш"""
    norm = get_delfia_norm(parameter, gestational_week, trimester)
    median = norm['median']
    
    if median > 0:
        mom = value / median
        
        if maternal_weight and parameter in ['PAPP_A', 'FREE_BETA_HCG', 'AFP', 'TOTAL_HCG']:
            weight_correction = np.sqrt(maternal_weight / 60)
            mom = mom / weight_correction
        
        return round(mom, 2)
    return 1.0

def get_age_multiplier(age, syndrome):
    """Ёшга кўра хавф кўпайтирувчисини олиш"""
    ages = sorted(AGE_MULTIPLIERS.keys())
    
    if age <= ages[0]:
        return AGE_MULTIPLIERS[ages[0]][syndrome]
    elif age >= ages[-1]:
        return AGE_MULTIPLIERS[ages[-1]][syndrome]
    else:
        for i in range(len(ages)-1):
            if ages[i] <= age <= ages[i+1]:
                low_age, high_age = ages[i], ages[i+1]
                low_mult = AGE_MULTIPLIERS[low_age][syndrome]
                high_mult = AGE_MULTIPLIERS[high_age][syndrome]
                
                # Линей интерполяция
                fraction = (age - low_age) / (high_age - low_age)
                return low_mult + fraction * (high_mult - low_mult)
    
    return 1.0

def calculate_syndrome_risks(age, nt_mom, papp_mom, hcg_mom, afp_mom=None, total_hcg_mom=None, ue3_mom=None):
    """Барча генетик синдромлар учун хавфларни ҳисоблаш"""
    
    risks = {}
    
    # Ёш хавфи
    age_risk_down = get_age_multiplier(age, 'downs')
    age_risk_edwards = get_age_multiplier(age, 'edwards')
    age_risk_patau = get_age_multiplier(age, 'patau')
    age_risk_turner = get_age_multiplier(age, 'turner')
    
    # Биринчи скрининг омиллари
    # Даун синдроми учун
    downs_risk = BASE_RISKS['downs'] * age_risk_down
    
    # PAPP-A омили
    if papp_mom < 0.3:
        downs_risk *= 3.0
    elif papp_mom < 0.4:
        downs_risk *= 2.0
    elif papp_mom < 0.5:
        downs_risk *= 1.5
    elif papp_mom > 2.5:
        downs_risk *= 1.2
    
    # Free β-hCG омили
    if hcg_mom < 0.2:
        downs_risk *= 2.5
    elif hcg_mom < 0.3:
        downs_risk *= 1.8
    elif hcg_mom > 2.5:
        downs_risk *= 2.0
    elif hcg_mom > 3.5:
        downs_risk *= 2.5
    
    # NT омили
    if nt_mom < 0.6:
        downs_risk *= 0.7
    elif nt_mom < 0.8:
        downs_risk *= 0.8
    elif nt_mom > 2.0:
        downs_risk *= 3.0
    elif nt_mom > 3.0:
        downs_risk *= 5.0
    elif nt_mom > 4.0:
        downs_risk *= 8.0
    
    risks['downs'] = min(downs_risk, 0.5)
    
    # Эдвардс синдроми учун
    edwards_risk = BASE_RISKS['edwards'] * age_risk_edwards
    
    # Эдвардс учун омиллар
    if papp_mom < 0.2:
        edwards_risk *= 4.0
    elif papp_mom < 0.3:
        edwards_risk *= 2.5
    
    if hcg_mom < 0.1:
        edwards_risk *= 3.0
    elif hcg_mom < 0.2:
        edwards_risk *= 2.0
    
    if nt_mom > 2.5:
        edwards_risk *= 4.0
    elif nt_mom > 3.0:
        edwards_risk *= 6.0
    
    risks['edwards'] = min(edwards_risk, 0.5)
    
    # Патау синдроми учун
    patau_risk = BASE_RISKS['patau'] * age_risk_patau
    
    if papp_mom < 0.2:
        patau_risk *= 5.0
    elif papp_mom < 0.3:
        patau_risk *= 3.0
    
    if hcg_mom < 0.15:
        patau_risk *= 3.5
    elif hcg_mom < 0.25:
        patau_risk *= 2.5
    
    if nt_mom > 2.8:
        patau_risk *= 5.0
    elif nt_mom > 3.5:
        patau_risk *= 8.0
    
    risks['patau'] = min(patau_risk, 0.5)
    
    # Тернер синдроми учун
    turner_risk = BASE_RISKS['turner'] * age_risk_turner
    
    if hcg_mom > 2.0:
        turner_risk *= 2.0
    elif hcg_mom > 3.0:
        turner_risk *= 3.0
    
    if nt_mom > 3.0:
        turner_risk *= 4.0
    elif nt_mom > 4.0:
        turner_risk *= 6.0
    
    risks['turner'] = min(turner_risk, 0.5)
    
    # Нейротубуляр дефект (НТД) учун
    ntd_risk = BASE_RISKS['ntd']
    
    if afp_mom and afp_mom > 2.5:
        ntd_risk = 0.01  # 1:100
    elif afp_mom and afp_mom > 2.0:
        ntd_risk = 0.02  # 1:50
    
    risks['ntd'] = min(ntd_risk, 0.5)
    
    # Ёш хавфи (алоҳида)
    risks['age_risk'] = {
        'downs': age_risk_down,
        'edwards': age_risk_edwards,
        'patau': age_risk_patau,
        'turner': age_risk_turner
    }
    
    # Иккиламчи скрининг омиллари (агар мавжуд бўлса)
    if all([afp_mom, total_hcg_mom, ue3_mom]):
        # Квад тест коррекцияси
        quad_correction = 1.0
        
        if afp_mom < 0.5:
            quad_correction *= 0.8
        elif afp_mom > 2.0:
            quad_correction *= 1.3
        
        if total_hcg_mom < 0.5:
            quad_correction *= 0.9
        elif total_hcg_mom > 2.0:
            quad_correction *= 1.8
        
        if ue3_mom < 0.5:
            quad_correction *= 1.5
        
        risks['downs'] *= quad_correction
        risks['edwards'] *= quad_correction * 1.2
        risks['patau'] *= quad_correction * 1.3
    
    return risks

def get_risk_category(risk_score):
    """Хавф категориясини аниқлаш"""
    if risk_score > 0.1:      # 1:10
        return "КРИТИК", "risk-critical", "#b71c1c"
    elif risk_score > 0.05:   # 1:20
        return "ЖУДА ЮҚОРИ", "risk-high", "#e65100"
    elif risk_score > 0.02:   # 1:50
        return "ЮҚОРИ", "risk-high", "#f57c00"
    elif risk_score > 0.01:   # 1:100
        return "ЎРТАЧА-ЮҚОРИ", "risk-medium", "#f57f17"
    elif risk_score > 0.005:  # 1:200
        return "ЎРТАЧА", "risk-medium", "#f9a825"
    elif risk_score > 0.001:  # 1:1000
        return "ПАСТ-ЎРТАЧА", "risk-low", "#388e3c"
    else:                     # 1:1000 дан кам
        return "ПАСТ", "risk-low", "#1b5e20"

def save_to_local_storage(data, filename):
    """Маълумотларни локал сақлаш"""
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        return True
    except Exception as e:
        st.error(f"Сақлашда хатолик: {e}")
        return False

def load_from_local_storage(filename):
    """Маълумотларни локалдан юклаш"""
    try:
        if os.path.exists(filename):
            with open(filename, 'r', encoding='utf-8') as f:
                return json.load(f)
        else:
            return []
    except Exception as e:
        print(f"Файлни ўқишда хатолик: {e}")
        return []

# ==================== КОНФИГУРАЦИЯ ====================
st.set_page_config(
    page_title="Генетик Синдромлар Хавф Бахолаш - DELFIA Revvity",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ==================== СТИЛЛАР ВА CSS ====================
st.markdown("""
<style>
    .main-header {
        font-size: 2.8rem;
        color: #0d47a1;
        text-align: center;
        margin-bottom: 1rem;
        font-weight: 800;
        background: linear-gradient(90deg, #0d47a1, #1565c0, #1976d2);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        padding: 25px;
        text-shadow: 3px 3px 8px rgba(13, 71, 161, 0.2);
        border-bottom: 5px solid #2196f3;
        border-radius: 12px;
        margin-top: 10px;
        border: 3px solid #bbdefb;
    }
    
    .sub-header {
        font-size: 1.6rem;
        color: #1565c0;
        text-align: center;
        margin-bottom: 2.5rem;
        font-weight: 600;
        background: linear-gradient(90deg, #e3f2fd, #bbdefb, #90caf9);
        padding: 20px;
        border-radius: 15px;
        border: 3px solid #2196f3;
        box-shadow: 0 8px 25px rgba(33, 150, 243, 0.2);
    }
    
    .syndrome-card {
        padding: 20px;
        border-radius: 15px;
        margin: 15px 0;
        border: 3px solid;
        box-shadow: 0 8px 20px rgba(0,0,0,0.1);
    }
    
    .downs-card { border-color: #ff6b6b; background: linear-gradient(135deg, #ffebee, #ffcdd2); }
    .edwards-card { border-color: #ff9800; background: linear-gradient(135deg, #fff3e0, #ffe0b2); }
    .patau-card { border-color: #ff5722; background: linear-gradient(135deg, #fbe9e7, #ffccbc); }
    .turner-card { border-color: #9c27b0; background: linear-gradient(135deg, #f3e5f5, #e1bee7); }
    .ntd-card { border-color: #4caf50; background: linear-gradient(135deg, #e8f5e9, #c8e6c9); }
    .age-risk-card { border-color: #2196f3; background: linear-gradient(135deg, #e3f2fd, #bbdefb); }
    
    .risk-critical {
        background: linear-gradient(135deg, #b71c1c, #d32f2f);
        color: white;
        padding: 15px 25px;
        border-radius: 25px;
        font-weight: bold;
        display: inline-block;
        border: 3px solid #ff5252;
        box-shadow: 0 6px 20px rgba(183, 28, 28, 0.3);
        animation: pulse 1.5s infinite;
        font-size: 1.2rem;
    }
    
    .risk-high {
        background: linear-gradient(135deg, #e65100, #f57c00);
        color: white;
        padding: 15px 25px;
        border-radius: 25px;
        font-weight: bold;
        display: inline-block;
        border: 3px solid #ffb74d;
        box-shadow: 0 6px 18px rgba(230, 81, 0, 0.3);
        font-size: 1.2rem;
    }
    
    .risk-medium {
        background: linear-gradient(135deg, #f57f17, #f9a825);
        color: #333;
        padding: 15px 25px;
        border-radius: 25px;
        font-weight: bold;
        display: inline-block;
        border: 3px solid #ffd54f;
        box-shadow: 0 6px 16px rgba(245, 127, 23, 0.3);
        font-size: 1.2rem;
    }
    
    .risk-low {
        background: linear-gradient(135deg, #1b5e20, #388e3c);
        color: white;
        padding: 15px 25px;
        border-radius: 25px;
        font-weight: bold;
        display: inline-block;
        border: 3px solid #66bb6a;
        box-shadow: 0 6px 16px rgba(27, 94, 32, 0.3);
        font-size: 1.2rem;
    }
    
    @keyframes pulse {
        0% { transform: scale(1); box-shadow: 0 0 0 0 rgba(183, 28, 28, 0.7); }
        50% { transform: scale(1.05); }
        70% { box-shadow: 0 0 0 15px rgba(183, 28, 28, 0); }
        100% { transform: scale(1); box-shadow: 0 0 0 0 rgba(183, 28, 28, 0); }
    }
    
    .metric-card {
        background: white;
        padding: 20px;
        border-radius: 15px;
        box-shadow: 0 5px 15px rgba(0,0,0,0.08);
        margin: 10px 0;
        border-left: 5px solid;
        transition: all 0.3s ease;
    }
    
    .metric-card:hover {
        transform: translateY(-3px);
        box-shadow: 0 8px 25px rgba(0,0,0,0.15);
    }
</style>
""", unsafe_allow_html=True)

# ==================== СЕССИЯ СОЗЛАМАЛАРИ ====================
if 'patient_id' not in st.session_state:
    st.session_state.patient_id = f"GEN-{datetime.now().strftime('%Y%m%d%H%M%S')}"
if 'patients_data' not in st.session_state:
    st.session_state.patients_data = load_from_local_storage("genetic_patients_data.json")
if 'current_patient' not in st.session_state:
    st.session_state.current_patient = {}
if 'screening_type' not in st.session_state:
    st.session_state.screening_type = "first"

# ==================== АСОСИЙ ИНТЕРФЕЙС ====================

# САРЛАВҲА
st.markdown('<h1 class="main-header">🧬 Генетик Синдромлар Хавф Бахолаш Дастури</h1>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Даун • Эдвардс • Патау • Тернер • НТД • Ёш хавфлари | DELFIA Revvity асосида</p>', unsafe_allow_html=True)

# СКРИНИНГ ТУРИ
st.markdown("### 📋 Скрининг турини танланг")
col1, col2 = st.columns(2)
with col1:
    if st.button("👶 Биринчи скрининг (10-14 ҳафта)", use_container_width=True):
        st.session_state.screening_type = "first"
        st.rerun()
with col2:
    if st.button("🤰 Иккиламчи скрининг (15-22 ҳафта)", use_container_width=True):
        st.session_state.screening_type = "second"
        st.rerun()

# САЙДБАР - БЕМОР МАЪЛУМОТЛАРИ
with st.sidebar:
    st.markdown("### 👤 Бемор маълумотлари")
    
    patient_name = st.text_input("Фамилия Исм Шариф", placeholder="Мадина Алиева")
    
    col_a, col_b = st.columns(2)
    with col_a:
        patient_age = st.number_input("Ёши", 15, 55, 30)
    with col_b:
        if st.session_state.screening_type == "first":
            gestational_age = st.number_input("Хомилалик (ҳафта)", 10, 14, 12)
        else:
            gestational_age = st.number_input("Хомилалик (ҳафта)", 15, 22, 18)
    
    height = st.number_input("Бўй (см)", 140, 200, 165)
    weight = st.number_input("Вазн (кг)", 40, 150, 65)
    
    if height > 0:
        bmi = calculate_bmi(weight, height)
        st.metric("📊 BMI", f"{bmi:.1f}")
    
    st.markdown("---")
    
    if st.session_state.screening_type == "first":
        st.markdown("### 🔬 Биринчи скрининг параметрлари")
        
        nt_measurement = st.slider("NT қалинлиги (мм)", 0.5, 10.0, 1.8, 0.1)
        
        papp_a_value = st.number_input(
            "PAPP-A Қиймати (U/L)", 
            0.1, 20.0, 1.4, 0.1
        )
        
        free_beta_hcg_value = st.number_input(
            "Free β-hCG Қиймати (ng/ml)", 
            1.0, 300.0, 80.0, 1.0
        )
    
    else:
        st.markdown("### 🔬 Иккиламчи скрининг параметрлари")
        
        afp_value = st.number_input(
            "AFP Қиймати (ng/ml)", 
            1.0, 200.0, 45.0, 1.0
        )
        
        total_hcg_value = st.number_input(
            "Total hCG Қиймати (IU/L)", 
            1000, 100000, 22000, 1000
        )
        
        ue3_value = st.number_input(
            "uE3 Қиймати (nmol/L)", 
            0.1, 20.0, 4.0, 0.1
        )
    
    st.markdown("---")
    calculate_btn = st.button("🧬 ГЕНЕТИК ХАВФЛАРНИ ҲИСОБЛАШ", 
                            type="primary", use_container_width=True)

# ==================== АСОСИЙ КОНТЕНТ ====================

if calculate_btn:
    if not patient_name:
        st.warning("⚠️ Илтимос, беморнинг исмини киритинг!")
    else:
        bmi = calculate_bmi(weight, height)
        
        with st.spinner("🧬 Генетик хавфлар ҳисобланади..."):
            if st.session_state.screening_type == "first":
                # Биринчи скрининг MoM ҳисоблаш
                papp_a_mom = calculate_mom_delfia(papp_a_value, 'PAPP_A', gestational_age, weight, "first")
                free_beta_hcg_mom = calculate_mom_delfia(free_beta_hcg_value, 'FREE_BETA_HCG', gestational_age, weight, "first")
                nt_mom = calculate_mom_delfia(nt_measurement, 'NT', gestational_age, weight, "first")
                
                # Хавфларни ҳисоблаш
                risks = calculate_syndrome_risks(
                    patient_age, nt_mom, papp_a_mom, free_beta_hcg_mom
                )
                
                # Маълумотларни сақлаш
                st.session_state.current_patient = {
                    'id': st.session_state.patient_id,
                    'name': patient_name,
                    'age': patient_age,
                    'screening_type': 'first',
                    'gestational_age': gestational_age,
                    'bmi': bmi,
                    'parameters': {
                        'nt': nt_measurement,
                        'nt_mom': nt_mom,
                        'papp_a': papp_a_value,
                        'papp_a_mom': papp_a_mom,
                        'free_beta_hcg': free_beta_hcg_value,
                        'free_beta_hcg_mom': free_beta_hcg_mom
                    },
                    'risks': risks,
                    'timestamp': datetime.now().isoformat()
                }
            
            else:
                # Иккиламчи скрининг MoM ҳисоблаш
                afp_mom = calculate_mom_delfia(afp_value, 'AFP', gestational_age, weight, "second")
                total_hcg_mom = calculate_mom_delfia(total_hcg_value, 'TOTAL_HCG', gestational_age, weight, "second")
                ue3_mom = calculate_mom_delfia(ue3_value, 'UE3', gestational_age, weight, "second")
                
                # Биринчи скрининг параметрларини сўраш (ихтиёрий)
                st.info("Биринчи скрининг параметрларини киритиш (ихтиёрий)")
                col_f1, col_f2, col_f3 = st.columns(3)
                
                with col_f1:
                    use_first_trimester = st.checkbox("Биринчи скрининг маълумотлари")
                
                if use_first_trimester:
                    first_gestational = st.number_input("Биринчи скрининг ҳафтаси", 10, 14, 12)
                    
                    col_p1, col_p2, col_p3 = st.columns(3)
                    with col_p1:
                        nt_measurement = st.number_input("NT (мм)", 0.5, 10.0, 1.8, 0.1)
                    with col_p2:
                        papp_a_value = st.number_input("PAPP-A (U/L)", 0.1, 20.0, 1.4, 0.1)
                    with col_p3:
                        free_beta_hcg_value = st.number_input("Free β-hCG (ng/ml)", 1.0, 300.0, 80.0, 1.0)
                    
                    papp_a_mom = calculate_mom_delfia(papp_a_value, 'PAPP_A', first_gestational, weight, "first")
                    free_beta_hcg_mom = calculate_mom_delfia(free_beta_hcg_value, 'FREE_BETA_HCG', first_gestational, weight, "first")
                    nt_mom = calculate_mom_delfia(nt_measurement, 'NT', first_gestational, weight, "first")
                    
                    # Хавфларни ҳисоблаш (икки скрининг билан)
                    risks = calculate_syndrome_risks(
                        patient_age, nt_mom, papp_a_mom, free_beta_hcg_mom,
                        afp_mom, total_hcg_mom, ue3_mom
                    )
                else:
                    # Фақат иккиламчи скрининг билан
                    risks = calculate_syndrome_risks(
                        patient_age, 1.0, 1.0, 1.0,
                        afp_mom, total_hcg_mom, ue3_mom
                    )
                
                # Маълумотларни сақлаш
                st.session_state.current_patient = {
                    'id': st.session_state.patient_id,
                    'name': patient_name,
                    'age': patient_age,
                    'screening_type': 'second',
                    'gestational_age': gestational_age,
                    'bmi': bmi,
                    'parameters': {
                        'afp': afp_value,
                        'afp_mom': afp_mom,
                        'total_hcg': total_hcg_value,
                        'total_hcg_mom': total_hcg_mom,
                        'ue3': ue3_value,
                        'ue3_mom': ue3_mom
                    },
                    'risks': risks,
                    'timestamp': datetime.now().isoformat()
                }
            
            st.session_state.patients_data.append(st.session_state.current_patient)
            save_to_local_storage(st.session_state.patients_data, "genetic_patients_data.json")
        
        st.success(f"✅ {patient_name} учун генетик хавфлар муваффақиятли ҳисобланди!")
        
        # БЕМОР МАЪЛУМОТЛАРИ
        st.markdown("### 📋 Бемор маълумотлари")
        col_info1, col_info2, col_info3, col_info4 = st.columns(4)
        
        with col_info1:
            st.metric("👤 Бемор", patient_name)
        with col_info2:
            st.metric("🎂 Ёши", f"{patient_age} йош")
        with col_info3:
            st.metric("🤰 Хомилалик", f"{gestational_age} ҳафта")
        with col_info4:
            st.metric("📊 BMI", f"{bmi:.1f}")
        
        # ГЕНЕТИК СИНДРОМЛАР ХАВФЛАРИ
        st.markdown("### 🧬 Генетик синдромлар хавфлари")
        
        # Даун синдроми
        with st.container():
            st.markdown('<div class="syndrome-card downs-card">', unsafe_allow_html=True)
            col_down1, col_down2, col_down3 = st.columns([2, 2, 3])
            
            with col_down1:
                st.markdown("#### Даун синдроми (Трисомия 21)")
                st.markdown("**Кифоялилик:** Интеллектуал нотўликлик, юрак аномалиялари")
            
            with col_down2:
                risk_value = risks['downs']
                risk_display = f"1:{int(1/risk_value)}" if risk_value > 0 else "1:∞"
                st.metric("Хавф нисбати", risk_display)
            
            with col_down3:
                category, risk_class, _ = get_risk_category(risk_value)
                st.markdown(f'<div class="{risk_class}">{category}</div>', unsafe_allow_html=True)
            
            st.markdown('</div>', unsafe_allow_html=True)
        
        # Эдвардс синдроми
        with st.container():
            st.markdown('<div class="syndrome-card edwards-card">', unsafe_allow_html=True)
            col_ed1, col_ed2, col_ed3 = st.columns([2, 2, 3])
            
            with col_ed1:
                st.markdown("#### Эдвардс синдроми (Трисомия 18)")
                st.markdown("**Кифоялилик:** Оғир кўп орган зарарланиш")
            
            with col_ed2:
                risk_value = risks['edwards']
                risk_display = f"1:{int(1/risk_value)}" if risk_value > 0 else "1:∞"
                st.metric("Хавф нисбати", risk_display)
            
            with col_ed3:
                category, risk_class, _ = get_risk_category(risk_value)
                st.markdown(f'<div class="{risk_class}">{category}</div>', unsafe_allow_html=True)
            
            st.markdown('</div>', unsafe_allow_html=True)
        
        # Патау синдроми
        with st.container():
            st.markdown('<div class="syndrome-card patau-card">', unsafe_allow_html=True)
            col_pa1, col_pa2, col_pa3 = st.columns([2, 2, 3])
            
            with col_pa1:
                st.markdown("#### Патау синдроми (Трисомия 13)")
                st.markdown("**Кифоялилик:** Ҳайвонот аномалиялари, НС зарарланиш")
            
            with col_pa2:
                risk_value = risks['patau']
                risk_display = f"1:{int(1/risk_value)}" if risk_value > 0 else "1:∞"
                st.metric("Хавф нисбати", risk_display)
            
            with col_pa3:
                category, risk_class, _ = get_risk_category(risk_value)
                st.markdown(f'<div class="{risk_class}">{category}</div>', unsafe_allow_html=True)
            
            st.markdown('</div>', unsafe_allow_html=True)
        
        # Тернер синдроми
        with st.container():
            st.markdown('<div class="syndrome-card turner-card">', unsafe_allow_html=True)
            col_tu1, col_tu2, col_tu3 = st.columns([2, 2, 3])
            
            with col_tu1:
                st.markdown("#### Тернер синдроми (45,X)")
                st.markdown("**Кифоялилик:** Бўй пастлиги, жинсий руксатсизлик")
            
            with col_tu2:
                risk_value = risks['turner']
                risk_display = f"1:{int(1/risk_value)}" if risk_value > 0 else "1:∞"
                st.metric("Хавф нисбати", risk_display)
            
            with col_tu3:
                category, risk_class, _ = get_risk_category(risk_value)
                st.markdown(f'<div class="{risk_class}">{category}</div>', unsafe_allow_html=True)
            
            st.markdown('</div>', unsafe_allow_html=True)
        
        # НТД (Нейротубуляр дефект)
        with st.container():
            st.markdown('<div class="syndrome-card ntd-card">', unsafe_allow_html=True)
            col_nt1, col_nt2, col_nt3 = st.columns([2, 2, 3])
            
            with col_nt1:
                st.markdown("#### Нейротубуляр дефект (НТД)")
                st.markdown("**Кифоялилик:** Спина бифида, анэнцефалия")
            
            with col_nt2:
                risk_value = risks['ntd']
                risk_display = f"1:{int(1/risk_value)}" if risk_value > 0 else "1:∞"
                st.metric("Хавф нисбати", risk_display)
            
            with col_nt3:
                category, risk_class, _ = get_risk_category(risk_value)
                st.markdown(f'<div class="{risk_class}">{category}</div>', unsafe_allow_html=True)
            
            st.markdown('</div>', unsafe_allow_html=True)
        
        # ЁШ ХАВФЛАРИ
        with st.container():
            st.markdown('<div class="syndrome-card age-risk-card">', unsafe_allow_html=True)
            st.markdown("#### 📊 Ёш бўйича хавф кўпайтирувчилари")
            
            age_risks = risks.get('age_risk', {})
            
            col_age1, col_age2, col_age3, col_age4 = st.columns(4)
            
            with col_age1:
                st.metric("Даун синдроми", f"{age_risks.get('downs', 1.0):.1f}x")
            with col_age2:
                st.metric("Эдвардс синдроми", f"{age_risks.get('edwards', 1.0):.1f}x")
            with col_age3:
                st.metric("Патау синдроми", f"{age_risks.get('patau', 1.0):.1f}x")
            with col_age4:
                st.metric("Тернер синдроми", f"{age_risks.get('turner', 1.0):.1f}x")
            
            st.markdown('</div>', unsafe_allow_html=True)
        
        # ТАҲЛИЛ ВА ГРАФИКЛАР
        st.markdown("### 📈 Хавф таҳлили")
        
        # Хавфлар тақсимоти
        col_chart1, col_chart2 = st.columns(2)
        
        with col_chart1:
            # Хавфлар диаграммаси
            syndromes = ['Даун', 'Эдвардс', 'Патау', 'Тернер', 'НТД']
            risk_values = [
                1/risks['downs'] if risks['downs'] > 0 else 10000,
                1/risks['edwards'] if risks['edwards'] > 0 else 10000,
                1/risks['patau'] if risks['patau'] > 0 else 10000,
                1/risks['turner'] if risks['turner'] > 0 else 10000,
                1/risks['ntd'] if risks['ntd'] > 0 else 10000
            ]
            
            fig = px.bar(
                x=syndromes,
                y=risk_values,
                title="Генетик синдромлар хавфлари (1:N)",
                labels={'x': 'Синдром', 'y': 'Хавф нисбати (1:N)'},
                color=syndromes,
                color_discrete_sequence=['#ff6b6b', '#ff9800', '#ff5722', '#9c27b0', '#4caf50']
            )
            fig.update_layout(height=400)
            st.plotly_chart(fig, use_container_width=True)
        
        with col_chart2:
            # Ёш хавфи графиги
            ages = list(AGE_MULTIPLIERS.keys())
            downs_mult = [AGE_MULTIPLIERS[age]['downs'] for age in ages]
            
            fig_age = px.line(
                x=ages,
                y=downs_mult,
                title="Ёш бўйича Даун синдроми хавфи",
                labels={'x': 'Ёш', 'y': 'Хавф кўпайтирувчиси'},
                markers=True
            )
            
            # Қўшимча синдромлар
            fig_age.add_scatter(
                x=ages,
                y=[AGE_MULTIPLIERS[age]['edwards'] for age in ages],
                mode='lines+markers',
                name='Эдвардс'
            )
            
            fig_age.add_scatter(
                x=ages,
                y=[AGE_MULTIPLIERS[age]['patau'] for age in ages],
                mode='lines+markers',
                name='Патау'
            )
            
            fig_age.update_layout(height=400, legend=dict(orientation="h", yanchor="bottom", y=1.02))
            st.plotly_chart(fig_age, use_container_width=True)
        
        # МАРКЕРЛАР ТАҲЛИЛИ
        st.markdown("### 🔬 Маркерлар таҳлили")
        
        if st.session_state.screening_type == "first":
            col_mark1, col_mark2, col_mark3 = st.columns(3)
            
            with col_mark1:
                st.metric("PAPP-A MoM", f"{papp_a_mom:.2f}")
                if papp_a_mom < 0.4:
                    st.error("Паст - хавф ошган")
                elif papp_a_mom > 2.5:
                    st.warning("Юқори - хавф ошган")
                else:
                    st.success("Нормал")
            
            with col_mark2:
                st.metric("Free β-hCG MoM", f"{free_beta_hcg_mom:.2f}")
                if free_beta_hcg_mom < 0.5:
                    st.error("Паст - хавф ошган")
                elif free_beta_hcg_mom > 2.0:
                    st.warning("Юқори - хавф ошган")
                else:
                    st.success("Нормал")
            
            with col_mark3:
                st.metric("NT MoM", f"{nt_mom:.2f}")
                if nt_measurement > 2.5:
                    st.error(f"Юқори: {nt_measurement} мм (норма: <2.5 мм)")
                else:
                    st.success(f"Нормал: {nt_measurement} мм")
        
        else:
            col_mark1, col_mark2, col_mark3 = st.columns(3)
            
            with col_mark1:
                st.metric("AFP MoM", f"{afp_mom:.2f}")
                if afp_mom < 0.5:
                    st.error("Паст - НТД хавфи")
                elif afp_mom > 2.0:
                    st.warning("Юқори - Даун хавфи")
                else:
                    st.success("Нормал")
            
            with col_mark2:
                st.metric("Total hCG MoM", f"{total_hcg_mom:.2f}")
                if total_hcg_mom < 0.5:
                    st.error("Паст - хавф ошган")
                elif total_hcg_mom > 2.0:
                    st.warning("Юқори - Даун хавфи")
                else:
                    st.success("Нормал")
            
            with col_mark3:
                st.metric("uE3 MoM", f"{ue3_mom:.2f}")
                if ue3_mom < 0.5:
                    st.error("Паст - Даун хавфи")
                else:
                    st.success("Нормал")
        
        # ТАВСИЯЛАР
        st.markdown("### 💡 Тиббий тавсиялар")
        
        with st.expander("#### 🏥 Хавф категориясига кўра тавсиялар", expanded=True):
            # Энг юқори хавфни аниқлаш
            max_risk = max(risks['downs'], risks['edwards'], risks['patau'], risks['turner'], risks['ntd'])
            max_syndrome = ""
            
            if max_risk == risks['downs']:
                max_syndrome = "Даун синдроми"
            elif max_risk == risks['edwards']:
                max_syndrome = "Эдвардс синдроми"
            elif max_risk == risks['patau']:
                max_syndrome = "Патау синдроми"
            elif max_risk == risks['turner']:
                max_syndrome = "Тернер синдроми"
            else:
                max_syndrome = "НТД"
            
            st.markdown(f"**Энг юқори хавф:** {max_syndrome} (1:{int(1/max_risk)})")
            
            if max_risk > 0.05:
                st.markdown("""
                **ШОШИЛИНЧ ЧОРАЛАР:**
                1. Дастурки генетик машварат (24 соат ичида)
                2. NIPT тести (но-инвазив пренатал тест)
                3. Амниоцентез ёки хорион биопсияси
                4. Фетал эхокардиография
                5. Ҳар ҳафта ультратовуш назорати
                """)
            elif max_risk > 0.01:
                st.markdown("""
                **ОЧИҚ ЧОРАЛАР:**
                1. Генетик машварат (72 соат ичида)
                2. Деталли ультратовуш таҳлили
                3. Қўшимча скрининг тестлари
                4. Ҳар 2 ҳафтада мониторинг
                """)
            elif max_risk > 0.001:
                st.markdown("""
                **НАЗОРАТ ЧОРАЛАРИ:**
                1. Генетик машварат (ихтиёрий)
                2. Мунтазам ультратовуш кўриқуви
                3. Парвардалик кўрсатмаларига риоя
                4. Ҳар 4-6 ҳафтада назорат
                """)
            else:
                st.markdown("""
                **НОРМАЛЬ ПАРВАРДАЛИК:**
                1. Стандарт скрининг дастури
                2. Регламент буйича ультратовуш
                3. Соглом турмуш тарзи
                4. Даво-профилактика витаминлари
                """)

else:
    st.markdown("""
    <div style="background: linear-gradient(135deg, #0d47a1 0%, #1976d2 100%); color: white; padding: 40px; border-radius: 20px; margin: 20px 0;">
        <h2 style="text-align: center; margin-bottom: 20px;">🧬 Генетик Синдромлар Хавф Бахолаш Дастурига Хуш Келибсиз!</h2>
        
        <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin-top: 30px;">
            <div style="background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px;">
                <h3>👶 Даун синдроми</h3>
                <p>Трисомия 21 - интеллектуал нотўликлик</p>
            </div>
            
            <div style="background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px;">
                <h3>⚠️ Эдвардс синдроми</h3>
                <p>Трисомия 18 - оғир кўп орган зарари</p>
            </div>
            
            <div style="background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px;">
                <h3>🔬 Патау синдроми</h3>
                <p>Трисомия 13 - неврологик аномалиялар</p>
            </div>
            
            <div style="background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px;">
                <h3>🧬 Тернер синдроми</h3>
                <p>45,X - жинсий хромосома аномалияси</p>
            </div>
            
            <div style="background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px;">
                <h3>📏 НТД</h3>
                <p>Нейротубуляр дефект - спина бифида</p>
            </div>
            
            <div style="background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px;">
                <h3>🎂 Ёш хавфи</h3>
                <p>Ёшга кўра хавф кўпайтирувчиси</p>
            </div>
        </div>
        
        <div style="text-align: center; margin-top: 40px;">
            <h3>📋 Дастурни ишлатиш учун:</h3>
            <p>1. Чеп томондаги панелда барча маълумотларни тўлдиринг</p>
            <p>2. Скрининг турини танланг (биринчи ёки иккиламчи)</p>
            <p>3. «ГЕНЕТИК ХАВФЛАРНИ ҲИСОБЛАШ» тугмасини босинг</p>
        </div>
    </div>
    """, unsafe_allow_html=True)

# ФУТЕР
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; padding: 20px;">
    <p style="font-size: 1.1rem; font-weight: bold; color: #0d47a1;">
        © 2024 Генетик Синдромлар Хавф Бахолаш Дастури | DELFIA Revvity асосида
    </p>
    <p style="font-size: 0.9rem; margin-top: 10px; color: #d32f2f;">
        ⚕️ ТИББИЙ ОГОҲЛАНТИРИШ: Бу дастур фақат ёрдамчи восита сифатида ишлатилади. 
        Ҳар қандай тиббий қарор қабул қилишдан олдин мутахассис шифокорга мурожаат қилинг.
    </p>
</div>
""", unsafe_allow_html=True)
