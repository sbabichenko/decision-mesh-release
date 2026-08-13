# The writeup's recommended configuration on Windows — PowerShell port
# of scripts/run_recommended.sh. Place decision_mesh.exe and the frozen
# data\smm_design_202606.csv next to this script (or pass paths).
#   Usage:  powershell -ExecutionPolicy Bypass -File run_recommended.ps1
param(
  [string]$Bin  = ".\decision_mesh.exe",
  [int]   $Seed = 7,
  [string]$Out  = ".\recommended_seed7",
  [string]$Data = ".\smm_design_202606.csv"
)
New-Item -ItemType Directory -Force -Path $Out | Out-Null

$env:DMESH_DATA = $Data
$env:DMESH_SPLIT = "1"
$env:DMESH_SPLIT_MODE = "local"
$env:DMESH_PL_SOLVER = "1"
$env:DMESH_HIER_FIT = "2"
$env:DMESH_HIER_MAP = "1"
$env:DMESH_QUAD_PRIOR = "1"
$env:DMESH_HEIGHT_LIMIT_LOGIT = "12"
$env:DMESH_IRLS_STEP_CAP = "4"
$env:DMESH_TAU_FLOOR_LOGIT_SD = "0.005"
$env:DMESH_STAGE_COUNT = "3"
$env:DMESH_MAX_GATE_ROUNDS = "80"
$env:DMESH_PARENT_VAR_SCALE = "0.25"
$env:DMESH_B3_USE_SCORE = "1"
$env:DMESH_B3_STAGE0_EXACT = "0"
# --- the recommended trio (writeup, Aug 5 2026 revision) ---
$env:DMESH_ONELAW = "0.1428,0.2415,0.0286,0.0038;8.5,8.1,5.4,0.0"
$env:DMESH_SIGMA0_MIN = "1.0"
$env:DMESH_SIGMA0_MAX = "6"
$env:DMESH_PI0_MIN = "0.05"
Remove-Item Env:DMESH_GATE_EXTRA_VAR, Env:DMESH_LAW_WEIGHTS, Env:DMESH_LAW_WEIGHTS_FINAL, Env:DMESH_CANDIDATE_BH, Env:DMESH_POINT_VAR -ErrorAction SilentlyContinue
$env:DMESH_DUMP = Join-Path $Out "run"

& $Bin 0 0.10 1 24 $Seed 1> (Join-Path $Out "stdout.log") 2> (Join-Path $Out "stderr.log")
Write-Host "done: $Out (held-out metrics: findstr HELDOUT $Out\stdout.log)"
