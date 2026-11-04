#!/usr/bin/env python3
"""Small OpenAI-compatible client used by GEMAgent review agents.

The runtime intentionally exposes only the subset required by the project:

    llm.bind(response_format={"type": "json_object"}).invoke(messages)

The returned object exposes ``.content`` as a string. Credentials are read from
``LLM_API_KEY``, ``LLM_BASE_URL``, and ``LLM_MODEL``.
"""

from __future__ import annotations

import asyncio
import os
import time

from dotenv import load_dotenv
from openai import OpenAI

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.dirname(os.path.dirname(_script_dir))
_workspace = os.path.dirname(_project_root)

for _env_dir in (
    _project_root,
    os.path.join(_workspace, "gemagent"),
    os.getcwd(),
    _script_dir,
):
    _env_file = os.path.join(_env_dir, ".env")
    if os.path.isfile(_env_file):
        load_dotenv(_env_file)
        break
else:
    load_dotenv()

api_key = os.getenv("LLM_API_KEY", "").strip()
base_url = os.getenv("LLM_BASE_URL", "").strip()
model_name = (
    os.getenv("LLM_MODEL_THINK", "").strip()
    or os.getenv("LLM_MODEL", "").strip()
    or "deepseek-chat"
)

_client: OpenAI | None = None
_LLM_USAGE_LOG = os.getenv("GEMAGENT_LLM_USAGE_LOG", "").strip()


def _get_client() -> OpenAI:
    """Create the shared API client on first use."""
    global _client
    if _client is not None:
        return _client
    if not api_key:
        raise RuntimeError(
            "LLM_API_KEY is not configured. Copy .env.example to .env and set "
            "the API credentials before running an LLM-backed review."
        )
    client_kwargs = {"api_key": api_key, "max_retries": 2, "timeout": 120}
    if base_url:
        client_kwargs["base_url"] = base_url
    _client = OpenAI(**client_kwargs)
    return _client


def _log_usage(line: str) -> None:
    print(line, flush=True)
    if _LLM_USAGE_LOG:
        with open(_LLM_USAGE_LOG, "a", encoding="utf-8") as handle:
            handle.write(line + "\n")


class _LlmResponse:
    """Minimal response object exposing ``.content``."""

    def __init__(self, content: str):
        self.content = content


def _convert_message(message) -> dict[str, str]:
    """Convert a dict or LangChain-style message into OpenAI chat format."""
    if isinstance(message, dict):
        return {
            "role": message.get("role", "user"),
            "content": message.get("content", ""),
        }

    role = getattr(message, "type", None)
    if role in ("human", "user"):
        role = "user"
    elif role in ("ai", "assistant"):
        role = "assistant"
    elif role == "system":
        role = "system"
    else:
        role = "user"

    content = getattr(message, "content", "")
    if not isinstance(content, str):
        content = str(content)
    return {"role": role, "content": content}


class _BoundLlm:
    """Bound model object exposing ``invoke(messages)``."""

    def __init__(
        self,
        model: str,
        response_format: dict | None = None,
        max_tokens: int = 16384,
        thinking: str = "enabled",
        max_empty_retries: int = 2,
    ):
        self._model = model
        self._response_format = response_format or {}
        self._max_tokens = max_tokens
        self._thinking = thinking
        self._max_empty_retries = int(max_empty_retries)

    def invoke(self, messages, **kwargs):
        del kwargs
        converted = [_convert_message(message) for message in messages]
        params = {
            "model": self._model,
            "messages": converted,
            "temperature": 0,
            "max_tokens": self._max_tokens,
            # The paper workflow enables reasoning. Callers may opt out with
            # bind(thinking="disabled").
            "extra_body": {"thinking": {"type": self._thinking}},
        }
        if self._response_format:
            params["response_format"] = self._response_format

        # JSON mode occasionally returns empty content. Retry a bounded number
        # of times before returning the empty response to the caller.
        attempts = self._max_empty_retries + 1
        for attempt in range(attempts):
            response = _get_client().chat.completions.create(**params)
            usage = getattr(response, "usage", None)
            if usage is not None:
                prompt_details = getattr(usage, "prompt_tokens_details", None)
                cached_tokens = getattr(prompt_details, "cached_tokens", 0) or 0
                _log_usage(
                    f"    [LLM] prompt={usage.prompt_tokens} "
                    f"cached={cached_tokens} completion={usage.completion_tokens}"
                )
            content = response.choices[0].message.content or ""
            if content.strip():
                return _LlmResponse(content)
            _log_usage(
                f"    [LLM] empty content, retry {attempt + 1}/{attempts}"
            )
            if attempt + 1 < attempts:
                time.sleep(1 + attempt)
        return _LlmResponse("")


class _Llm:
    """Lightweight ChatOpenAI-compatible ``bind().invoke()`` facade."""

    def __init__(self, model: str):
        self._model = model
        self.model_name = model
        # Keep LangChain-style attributes for diagnostic scripts.
        self.openai_api_base = base_url
        self.openai_api_key = api_key
        self.temperature = 0
        self.max_tokens = 16384
        self.max_retries = 3

    def bind(self, **kwargs):
        return _BoundLlm(
            self._model,
            kwargs.get("response_format"),
            kwargs.get("max_tokens", 16384),
            kwargs.get("thinking", "enabled"),
            kwargs.get("max_empty_retries", 2),
        )

    def invoke(self, messages, **kwargs):
        return self.bind(**kwargs).invoke(messages)

    async def ainvoke(self, messages, **kwargs):
        """Async facade used by review steps that call ``await llm.ainvoke()``."""
        return await asyncio.to_thread(self.invoke, messages, **kwargs)


llm = _Llm(model_name)
