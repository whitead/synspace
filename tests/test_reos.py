from synspace.reos import REOS


def test_can_fetch_file():
    """Test that the REOS class can fetch the alert_collection.csv file."""
    reos = REOS()
    assert reos.rule_path is not None
