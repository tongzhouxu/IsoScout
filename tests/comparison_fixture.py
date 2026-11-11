"""Minimal assembly runner result for orchestration tests."""
def comparison_result(policy=None):
    return {"commands": [], "comparisons": [], "alignment_jobs": 1, "cache_hits": 0,
            "elapsed_seconds": 0., "status": "PASS", "tools": {}, "policy": policy or {},
            "policy_sha256": "fixture", "implementation_sha256": "fixture", "cache_directory": "fixture"}
