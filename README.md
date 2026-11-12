# VMEC Websites

VMEC-lyh 程序手册和 NCSX 自由边界算例的公开网站，无需登录。

| 目录 | 内容 | 网站 |
| --- | --- | --- |
| [Source-Code](Source-Code/README.md) | 程序手册、改动说明、运行参数、验证记录和 [VMEC-lyh 源码快照](Source-Code/VMEC-lyh/) | [程序手册](https://yihui-l.github.io/VMEC-websites/Source-Code/) |
| [VMEC-lyh-example](VMEC-lyh-example/README.md) | NCSX 算例的输入、两版 WOUT、数据对比和一致性检查 | [算例网站](https://yihui-l.github.io/VMEC-websites/VMEC-lyh-example/) |

```text
VMEC-websites/
├── Source-Code/
│   ├── VMEC-lyh/
│   └── 程序手册、版本说明与网页资源
├── VMEC-lyh-example/
│   └── NCSX/
├── README.md
└── .nojekyll
```

程序手册迁自原 `Documents/`。源码快照的版本及许可见 [源码版本说明](Source-Code/SOURCE_VERSION.md)。NCSX 资料保留原始输入和 WOUT，计算条件、诊断定义和文件校验和见 [数据对比记录](VMEC-lyh-example/NCSX/wout_comparison.md)。

GitHub Pages 从 `gh-pages` 分支的根目录发布。更新网站时，将手册和算例的改动同步到该分支。两个站点入口分别位于 `Source-Code/` 和 `VMEC-lyh-example/`。手册中的源码链接指向本仓库的完整源码快照。
