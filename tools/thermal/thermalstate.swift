import Foundation
// Prints macOS's own thermal pressure assessment: 0 nominal, 1 fair, 2 serious, 3 critical
print(ProcessInfo.processInfo.thermalState.rawValue)
