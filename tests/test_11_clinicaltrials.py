"""
tests/test_11_clinicaltrials.py
Reproduce bugs reportados en 11_clinicaltrials_pubmed.py (plan B3).
Ambiente: omics-py
Ejecutar: python tests/test_11_clinicaltrials.py
"""
import sys
import re
import importlib.util
from pathlib import Path

# Agregar scripts/ al path para importar desde el script
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

# ─── Importar solo las funciones/constantes necesarias (sin ejecutar main) ───
# Cargar el módulo sin ejecutar __main__
spec = importlib.util.spec_from_file_location(
    "script11",
    Path(__file__).resolve().parent.parent / "scripts" / "11_clinicaltrials_pubmed.py"
)
mod = importlib.util.module_from_spec(spec)
# Parchear __name__ para evitar que corra main()
mod.__name__ = "script11"
spec.loader.exec_module(mod)

PASS = []
FAIL = []

def check(name, condition, msg_fail=""):
    if condition:
        PASS.append(name)
        print(f"  PASS  {name}")
    else:
        FAIL.append(name)
        print(f"  FAIL  {name}  —  {msg_fail}")

# =============================================================================
# Bug 1: HNSCC_PUBMED_TERM contiene "squamous cell carcinoma"[Title/Abstract]
#        sin calificador "head and neck" → captura SCC de otros sitios
# =============================================================================
print("\n--- Bug 1: HNSCC_PUBMED_TERM sobreconteo ---")
term = mod.HNSCC_PUBMED_TERM
# El term NO debe incluir "squamous cell carcinoma" sin modificador de head/neck
# Comprobación: buscar el patrón problemático
scc_bare = re.search(r'"squamous cell carcinoma"\[Title/Abstract\]', term)
check(
    "pubmed_term_no_bare_SCC",
    scc_bare is None,
    f'HNSCC_PUBMED_TERM contiene "squamous cell carcinoma"[Title/Abstract] sin '
    f'calificador head/neck — capturará SCC de pulmón/cérvix/esófago. '
    f'Valor actual: {term}'
)

# El term SÍ debe exigir "head and neck" de alguna forma
has_hn = ("head and neck" in term.lower() or "head neck" in term.lower()
          or "HNSCC" in term or "hnscc" in term.lower())
check(
    "pubmed_term_requires_head_neck",
    has_hn,
    f'HNSCC_PUBMED_TERM no contiene calificador head/neck. Valor: {term}'
)

# =============================================================================
# Bug 2: is_hnscc usa substring crudo — "oral cancer" como kw puede matchear
#        condiciones no-HNSCC que contengan "oral" (ej: "transoral approach")
# =============================================================================
print("\n--- Bug 2: is_hnscc substring crudo ---")

# Simular la lógica actual: any(kw.lower() in conditions_str for kw in HNSCC_CT_TERMS)
def is_hnscc_current(conditions):
    conditions_str = " ".join(c.lower() for c in conditions)
    return any(kw.lower() in conditions_str for kw in mod.HNSCC_CT_TERMS)

# Caso problemático: "transoral robotic surgery complications" no es HNSCC como condición
# pero "oral" está en la string de condiciones
false_pos_case = ["Transoral approach complications", "Postoperative complications"]
result_fp = is_hnscc_current(false_pos_case)
check(
    "is_hnscc_no_false_positive_transoral",
    not result_fp,
    f'is_hnscc retorna True para {false_pos_case} por substring "oral" — falso positivo'
)

# Caso verdadero positivo: debe seguir funcionando
true_pos_case = ["Head and Neck Squamous Cell Carcinoma", "Oropharyngeal Cancer"]
result_tp = is_hnscc_current(true_pos_case)
check(
    "is_hnscc_true_positive_hnscc",
    result_tp,
    f'is_hnscc retorna False para caso HNSCC real — falso negativo'
)

# =============================================================================
# Bug 3: drug_confirmed = True cuando interventions está vacío
# =============================================================================
print("\n--- Bug 3: drug_confirmed infla conteos ---")

# Verificar que el código fuente del script NO contiene el patrón buggy
import inspect
src = inspect.getsource(mod.query_clinicaltrials)

buggy_pattern = "if interventions else True"
check(
    "drug_confirmed_false_when_no_interventions",
    buggy_pattern not in src,
    f'query_clinicaltrials aún contiene el patrón "{buggy_pattern}" — '
    f'asume drug_confirmed=True cuando interventions está vacío'
)

# Verificar que el fix sí está (rechazar cuando no hay intervenciones)
fix_pattern = "drug_confirmed = False"
check(
    "drug_confirmed_fix_present",
    fix_pattern in src,
    f'Fix no encontrado en query_clinicaltrials: espera "{fix_pattern}" cuando interventions vacío'
)

# =============================================================================
# Bug 4: Input file apunta a 10_top20_candidates.tsv (no existe en pipeline actual)
# =============================================================================
print("\n--- Bug 4: Input file path obsoleto ---")

# Después del fix, el script usa INPUT_SCORED en vez de INPUT_TOP20
has_scored = hasattr(mod, "INPUT_SCORED")
check(
    "input_uses_scored_not_top20",
    has_scored and not hasattr(mod, "INPUT_TOP20"),
    f'Script aún usa INPUT_TOP20 (10_top20_candidates.tsv) en vez de INPUT_SCORED '
    f'(10_all_candidates_scored.tsv + LOD filter)'
)

# El archivo al que apunta debe existir
if has_scored:
    check(
        "input_file_exists",
        mod.INPUT_SCORED.exists(),
        f'Archivo input no encontrado: {mod.INPUT_SCORED}'
    )

# =============================================================================
# Resumen
# =============================================================================
print(f"\n{'='*50}")
print(f"PASS: {len(PASS)}  FAIL: {len(FAIL)}")
if FAIL:
    print(f"Bugs confirmados: {', '.join(FAIL)}")
    sys.exit(1)
else:
    print("TEST PASSED — todos los bugs corregidos")
    sys.exit(0)
