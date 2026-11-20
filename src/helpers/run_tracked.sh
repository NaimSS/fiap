# runs a solver with a timeout and a memory cap. source it from a bench script:
#   source "$ROOT/src/helpers/run_tracked.sh"
#   run_tracked <timeout_s> <mem_cap_mb> <cmd> [args...] < input
#
# timeout: perl alarm + exec, the solver itself gets SIGALRM (exit 142).
# memory: MEM_CAP_MB, only works on binaries built with mem_cap.h.
#
# sets RT_EXIT, RT_STATUS (OK | TIMEOUT | OOM | ASSERT | ERROR(<rc>)),
# RT_MEM_MB (peak from mem_cap.h), RT_RSS_MB (max RSS from /usr/bin/time),
# RT_OUT and RT_ERR (files with stdout / stderr; rt_cleanup removes them).

RT_OUT="$(mktemp -t rt_out.XXXXXX)"
RT_ERR="$(mktemp -t rt_err.XXXXXX)"
rt_cleanup(){ rm -f "$RT_OUT" "$RT_ERR"; }

# set before the first call, for callers with set -u
RT_EXIT=0; RT_STATUS=""; RT_MEM_MB=""; RT_RSS_MB=""

# max-RSS reporter: BSD time -l (macOS, bytes) or GNU time -v (Linux, KiB)
_RT_TIME=(); _RT_RSS_UNIT=0
if [[ -x /usr/bin/time ]]; then
    if /usr/bin/time -l true >/dev/null 2>&1; then
        _RT_TIME=(/usr/bin/time -l); _RT_RSS_UNIT=1
    elif /usr/bin/time -v true >/dev/null 2>&1; then
        _RT_TIME=(/usr/bin/time -v); _RT_RSS_UNIT=1024
    fi
fi

run_tracked(){
    local tmo="$1" cap_mb="$2"; shift 2
    RT_EXIT=0; RT_STATUS=""; RT_MEM_MB=""; RT_RSS_MB=""
    # alarm(0) would silently disable the timeout
    if ! [[ "$tmo" =~ ^[0-9]+$ ]] || (( tmo < 1 )); then
        echo "run_tracked: timeout must be an integer >= 1 s (got '$tmo')" >&2
        RT_EXIT=2; RT_STATUS="ERROR(args)"; : > "$RT_OUT"; : > "$RT_ERR"
        return 0
    fi

    # 2>/dev/null hides bash's "Alarm clock" message; the `exit $?` keeps bash
    # from exec-ing the command directly, which would bring it back
    ( MEM_CAP_MB="$cap_mb" MEM_CAP_REPORT=1 \
        ${_RT_TIME[@]+"${_RT_TIME[@]}"} perl -e 'alarm shift @ARGV; exec @ARGV or exit 127' \
        "$tmo" "$@" > "$RT_OUT" 2> "$RT_ERR"; exit $? ) 2>/dev/null
    RT_EXIT=$?

    if   (( RT_EXIT == 0 ));   then RT_STATUS=OK
    elif (( RT_EXIT == 142 )); then RT_STATUS=TIMEOUT
    elif grep -q "memory cap exceeded" "$RT_ERR"; then RT_STATUS=OOM
    elif (( RT_EXIT == 134 )); then RT_STATUS=ASSERT
    else RT_STATUS="ERROR($RT_EXIT)"
    fi

    local kb
    kb=$(grep -o 'mem_cap_peak_kb: [0-9]*' "$RT_ERR" | tail -1 | awk '{print $2}')
    [[ -n "$kb" ]] && RT_MEM_MB=$(awk -v k="$kb" 'BEGIN{printf "%.3f", k/1024}')

    local rss
    if (( _RT_RSS_UNIT == 1 )); then       # BSD: "<bytes>  maximum resident set size"
        rss=$(awk '/maximum resident set size/{print $1}' "$RT_ERR" | tail -1)
    else                                   # GNU: "Maximum resident set size (kbytes): <KiB>"
        rss=$(awk -F': *' '/Maximum resident set size/{print $2}' "$RT_ERR" | tail -1)
    fi
    [[ -n "$rss" ]] && RT_RSS_MB=$(awk -v b="$rss" -v u="$_RT_RSS_UNIT" 'BEGIN{printf "%.1f", b*u/1048576}')
    return 0
}

# reads the result line of the tree solvers from $RT_OUT into TD_VALUE,
# TD_TIME_S, TD_CUTS, TD_BS_STEPS, TD_SOLVER_STATUS. an OK run becomes
# ERROR(parse) if something is missing, or GAP_OR_COLLINEAR if not exact
rt_parse_tree_dp(){
    TD_VALUE=""; TD_TIME_S=""; TD_CUTS=""; TD_BS_STEPS=""; TD_SOLVER_STATUS=""
    local line
    line=$(grep -m1 -- '  - k=' "$RT_OUT")
    TD_VALUE=$(awk '{print $1}' <<< "$line")
    TD_TIME_S=$(sed -nE 's/.*, time: ([0-9.eE+-]+)s.*/\1/p' <<< "$line")
    TD_CUTS=$(sed -nE 's/.*, cuts: ([0-9]+).*/\1/p' <<< "$line")
    TD_BS_STEPS=$(sed -nE 's/.*, bs_steps: ([0-9]+).*/\1/p' <<< "$line")
    TD_SOLVER_STATUS=$(sed -nE 's/.*, status: ([A-Za-z_-]+).*/\1/p' <<< "$line")
    [[ "$RT_STATUS" == OK ]] || return 0
    if [[ -z "$TD_VALUE" || -z "$TD_TIME_S" || -z "$TD_SOLVER_STATUS" ]]; then
        RT_STATUS="ERROR(parse)"
    elif [[ "$TD_SOLVER_STATUS" != exact* ]]; then
        RT_STATUS="$TD_SOLVER_STATUS"
    fi
}
