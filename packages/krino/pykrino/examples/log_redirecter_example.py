#!/usr/bin/env python3
"""
Example script demonstrating LogRedirecter functionality in pykrino.

This script shows how to:
1. Create a LogRedirecter to capture krino log output
2. Use the captured logs for debugging and testing
3. Clear the log buffer when needed

Note: This is a conceptual example. Actual usage requires:
- Krino operations that generate log output
- Proper error handling
"""

import sys
import os

# Add pykrino to path if not installed
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    import pykrino
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure pykrino is built.")
    sys.exit(1)


def example_basic_log_redirecter():
    """Example of basic LogRedirecter usage"""
    print("=== Basic LogRedirecter Example ===")

    # Create a LogRedirecter
    # This will capture all krino log output from this point forward
    redirecter = pykrino.LogRedirecter()
    print("Created LogRedirecter")

    # Get current log output (should be empty initially)
    log_output = redirecter.get_log()
    print(f"Current log output length: {len(log_output)} characters")

    # Clear the log buffer
    redirecter.clear()
    print("Cleared log buffer")

    # Get log output after clearing
    log_output = redirecter.get_log()
    print(f"Log output after clear: {len(log_output)} characters")


def example_log_capture():
    """Example of capturing logs during krino operations"""
    print("\n=== Log Capture Example ===")

    # Create a LogRedirecter
    redirecter = pykrino.LogRedirecter()
    print("Created LogRedirecter for capturing logs")

    # At this point, any krino operations would have their output
    # captured by the redirecter. For example:
    #
    # mesh = pykrino.MeshFromFile("example_mesh.e")
    # mesh.populate_mesh()
    #
    # Any log messages from these operations would be captured

    # Get the captured logs
    captured_logs = redirecter.get_log()
    print(f"Captured {len(captured_logs)} characters of log output")

    # Clear the buffer for next operation
    redirecter.clear()


def example_log_analysis():
    """Example of analyzing captured logs"""
    print("\n=== Log Analysis Example ===")

    # Create a LogRedirecter
    redirecter = pykrino.LogRedirecter()

    # After performing krino operations, analyze the logs
    logs = redirecter.get_log()

    # Example analysis:
    if logs:
        print(f"Log contains {len(logs)} characters")
        print(f"Log has {len(logs.splitlines())} lines")

        # Check for specific patterns
        if "warning" in logs.lower():
            print("Warning detected in logs")
        if "error" in logs.lower():
            print("Error detected in logs")
    else:
        print("No log output captured")

    redirecter.clear()


def example_multiple_operations():
    """Example of capturing logs from multiple operations"""
    print("\n=== Multiple Operations Example ===")

    # Create a LogRedirecter
    redirecter = pykrino.LogRedirecter()

    # Operation 1
    # mesh1 = pykrino.MeshFromFile("mesh1.e")
    # mesh1.populate_mesh()
    logs_after_op1 = redirecter.get_log()
    print(f"After operation 1: {len(logs_after_op1)} characters")

    # Clear for next operation
    redirecter.clear()

    # Operation 2
    # mesh2 = pykrino.MeshFromFile("mesh2.e")
    # mesh2.populate_mesh()
    logs_after_op2 = redirecter.get_log()
    print(f"After operation 2: {len(logs_after_op2)} characters")

    # Compare logs
    print(f"Operation 1 generated {len(logs_after_op1)} chars")
    print(f"Operation 2 generated {len(logs_after_op2)} chars")


def main():
    """Main function to run examples"""
    print("PyKrino LogRedirecter Examples")
    print("=" * 60)

    # Run examples
    example_basic_log_redirecter()
    example_log_capture()
    example_log_analysis()
    example_multiple_operations()

    print("\n" + "=" * 60)
    print("All LogRedirecter examples completed!")
    print("\nKey Takeaways:")
    print("1. LogRedirecter captures krino log output to a buffer")
    print("2. Use get_log() to retrieve captured output as a string")
    print("3. Use clear() to reset the log buffer")
    print("4. Logs are automatically captured once LogRedirecter is created")
    print("5. LogRedirecter is useful for testing and debugging")


if __name__ == "__main__":
    main()
