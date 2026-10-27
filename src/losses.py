import torch
import torch.nn as nn
import torch.nn.functional as F


class PositionWeightedCrossEntropy(nn.Module):
    def __init__(self, position_class_weights):
        super().__init__()
        self.register_buffer(
            "position_class_weights",
            position_class_weights.float(),
        )

    def forward(self, logits, targets):
        targets = targets.long()

        per_position_loss = F.cross_entropy(
            logits,
            targets,
            reduction="none",
        )

        weights = self.position_class_weights.T
        target_weights = weights[
            targets,
            torch.arange(9, device=targets.device),
        ]

        return (per_position_loss * target_weights).sum() / target_weights.sum().clamp_min(1e-12)